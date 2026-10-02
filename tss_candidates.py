#!/usr/bin/env python

"""
Write every TSS candidate (all enriched 5' peaks, both strands) to one table,
without applying the TSS threshold.

Candidates with a control peak within +/- the peak alignment error get the R9
enrichment ratio (RPM enriched / RPM control). Candidates without one are kept
and flagged "enriched_only"; their ratio uses the raw control 5' read count
around the TSS, floored at a pseudocount, so every candidate has a finite ratio.
"""

import argparse
import numpy as np
import pandas as pd

STRANDS = {"plus": "+", "minus": "-"}

COLUMNS = ["candidate_id", "chrom", "strand", "tss_pos", "peak_start", "peak_end",
	"enriched_reads", "enriched_rpm", "control_match", "control_pos", "control_reads",
	"control_rpm", "control_window_reads", "enrichment_ratio", "enrichment_ratio_all",
	"above_threshold"]


"""
Read a clustered peak file written by peak_clustering.r
"""
def readClustered(file):
	df = pd.read_csv(file, index_col=0)
	df["chr"] = df["chr"].astype(str)
	return df


"""
Sum the 5' read ends of a bedgraph over 1-based positions [pos - error, pos + error]
"""
class WindowCounter:
	def __init__(self, bedgraph):
		df = pd.read_csv(bedgraph, sep='\t', header=None, usecols=[0, 1, 2, 3],
			names=["chrom", "start", "end", "count"], dtype={"chrom": str})
		self.by_chrom = {}
		for chrom, d in df.groupby("chrom"):
			d = d.sort_values("start")
			starts = d["start"].to_numpy()
			ends = d["end"].to_numpy()
			counts = d["count"].to_numpy()
			cum = np.concatenate([[0], np.cumsum(counts * (ends - starts))[:-1]])
			self.by_chrom[chrom] = (starts, ends, counts, cum)

	# number of read ends on 0-based bases < x
	def _cumulative(self, chrom, x):
		starts, ends, counts, cum = self.by_chrom[chrom]
		i = np.searchsorted(starts, x, side="right") - 1
		out = np.zeros(len(x))
		ok = i >= 0
		j = i[ok]
		out[ok] = cum[j] + counts[j] * (np.minimum(x[ok], ends[j]) - starts[j])
		return out

	def count(self, chrom, pos, error):
		if chrom not in self.by_chrom:
			return np.zeros(len(pos))
		pos = np.asarray(pos)
		return self._cumulative(chrom, pos + error) - self._cumulative(chrom, pos - error - 1)


def candidatesForStrand(enriched, control, bedgraph, sign, error, pseudocount, threshold):
	symb = STRANDS[sign]
	enr = readClustered(enriched)
	ctl = readClustered(control)
	counter = WindowCounter(bedgraph)

	if ctl.empty:
		control_mapped = np.nan
	else:
		control_mapped = ctl["mapped_reads"].iloc[0]

	rows = []
	for _, e in enr.iterrows():
		pos = int(e["decision_v"])
		same = ctl[(ctl["chr"] == e["chr"]) & ((ctl["decision_v"] - pos).abs() <= error)]
		row = {
			"candidate_id": "{}_{}_{}".format(e["chr"], sign, pos),
			"chrom": e["chr"],
			"strand": symb,
			"tss_pos": pos,
			"peak_start": int(e["start_peak"]),
			"peak_end": int(e["end_peak"]),
			"enriched_reads": e["full_cov_clust"],
			"enriched_rpm": e["RPM"],
		}
		if same.empty:
			row["control_match"] = "enriched_only"
		else:
			# nearest control peak; on a tie the one with most reads
			same = same.assign(dist=(same["decision_v"] - pos).abs()) \
				.sort_values(["dist", "full_cov_clust"], ascending=[True, False])
			c = same.iloc[0]
			row.update({
				"control_match": "matched",
				"control_pos": int(c["decision_v"]),
				"control_reads": c["full_cov_clust"],
				"control_rpm": c["RPM"],
				"enrichment_ratio": e["RPM"] / c["RPM"],
			})
		rows.append(row)

	df = pd.DataFrame(rows, columns=COLUMNS)
	if df.empty:
		return df

	for chrom, idx in df.groupby("chrom").groups.items():
		df.loc[idx, "control_window_reads"] = counter.count(chrom, df.loc[idx, "tss_pos"].to_numpy(), error)

	for col in ["control_pos", "control_reads", "control_window_reads"]:
		df[col] = df[col].round().astype("Int64")

	floor_rpm = 1000000 * np.maximum(df["control_window_reads"].astype(float), pseudocount) / control_mapped
	df["enrichment_ratio_all"] = df["enrichment_ratio"].where(df["control_match"] == "matched",
		df["enriched_rpm"] / floor_rpm)
	df["above_threshold"] = (df["control_match"] == "matched") & (df["enrichment_ratio"] > threshold)
	return df


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description=__doc__)
	for sign in STRANDS:
		parser.add_argument("--enriched-" + sign, required=True)
		parser.add_argument("--control-" + sign, required=True)
		parser.add_argument("--control-bedgraph-" + sign, required=True)
	parser.add_argument("--error", type=int, required=True, help="peak alignment error (nt)")
	parser.add_argument("--pseudocount", type=float, default=1, help="control read floor for enriched-only peaks")
	parser.add_argument("--threshold", type=float, required=True, help="TSS Threshold, only used for the above_threshold flag")
	parser.add_argument("--output", required=True)
	args = vars(parser.parse_args())

	tables = [candidatesForStrand(args["enriched_" + sign], args["control_" + sign], args["control_bedgraph_" + sign],
		sign, args["error"], args["pseudocount"], args["threshold"]) for sign in STRANDS]
	out = pd.concat([t for t in tables if not t.empty], ignore_index=True) if any(not t.empty for t in tables) \
		else pd.DataFrame(columns=COLUMNS)
	out = out.sort_values(["chrom", "tss_pos", "strand"])
	out.to_csv(args["output"], sep='\t', index=False)

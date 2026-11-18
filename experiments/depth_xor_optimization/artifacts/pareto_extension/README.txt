Pareto comparison against the existing xordepth evidence for four fixed matrices.

No algorithm source changes. Only independent synthesis, native verification, read-only exact audits and plotting.
V1: OUR_DCCSE,16 greedy restarts,seed1.
V2: OUR_DCCSE_V2,8 restarts per strategy (16 total),seed1,default refinements.
Neither loads --incumbent or any repository circuit.
Every new run has an external180s wall-clock cap. All completed successfully.

New caps: Vandermonde/Cauchy7,8; GF8 8,9; RS7,8,9,11,14.
Reused caps: Vandermonde/RS/Cauchy5,6; GF8 6,7, from the previous verified V2 experiment. There are19 completed circuits per independent method (38 total), including22 new runs.

Repository group includes all60 candidates of the original four frontier experiments,20 entries from freshly rerun native --depth-aware --share2,3,0 frontiers, and4 previously completed verified native bounded results. Historical direct CSE, MW_SQ, UBER_XSET, legacy Uber and ring pipelines are imported exactly through the same xordepth --also evidence used earlier; original and rebalanced variants remain separate rows. No historical solver is silently described as a fresh equal-time run. Failed/timeout bounded records have no circuit point and are not treated as infeasibility proofs.
This covers the existing four-matrix evidence and supported variants above, not every possible parameterization or a new exhaustive search of every backend. Historical search budgets differ from our runs.

All122 circuit records are independently checked for exact GF2 target outputs, valid dependencies, live-gate counts and actual maximum output depth. Our38 results also have successful native verification. The20 fresh depth-aware candidates share4 successful native frontier verifications. Some repository records duplicate circuits or metrics across methods; Pareto extraction deduplicates metric points.

Pareto definition: (gates,actual depth) is retained iff no other candidate has both no more gates and no greater depth, with at least one strict improvement. Selection uses actual depth, not the requested cap. A candidate discovered under a looser cap can be reused at any cap covering its actual depth. All selected circuits can be traced to all_candidates.csv and source JSON evidence.

Files:
summary.json: per-method and combined observed Pareto fronts.
all_candidates.csv: every original and added circuit record (including dominated ones),122 rows.
best_under_depth.csv: best observed feasible result under each displayed integer cap, reusing all completed circuits. RS10/12/13 were not fresh synthesis caps; these entries select known feasible circuits.
report.html: Chinese report, all candidate rows grouped by matrix, graph and Pareto tables.
pareto_comparison.png/.pdf: scientific plots; steps are best-known feasible costs, not proof of optimality or interpolation claims.
reused_evidence/: prior successful circuits, inputs and original frontier evidence.
run*.sh: native synthesis and depth-aware reproduction commands; set REPO_ROOT to the already patched checkout.
analyze.py: only verification/reporting, no synthesis implementation. It can use sibling four_matrices/ and v2_results/ from the previous artifact, or fall back to its bundled reused_evidence/ folder.

Main finding: no additional lower-gate point appeared at the expanded caps with these fixed strategies/seeds. V2 dominates listed repository fronts for Cauchy and GF8; Vandermonde still has the repository203/6 tradeoff; RS has our145/5 shallow option and much smaller repository circuits at deeper depths. These are observed finite candidate fronts, not global optimal Pareto fronts.

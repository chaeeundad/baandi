Architect decisions (rule R2 stays frozen; these only clarify it):
1. si32: use the ticket's special count, 5 standard starts and 1 strict reference per variant (cost). All other systems: 10 standard + 2 strict.
2. Total cost: judge by measured single-thread wall time of the whole SCF (median over starts, same machine and build, paired free_electron vs atomic, alternating order), which includes the atomic start honestly. Report H applications and MAC counts as diagnostics. The +5 percent limit applies to median wall time. Mac numbers are observations; dev08 confirms later.
Proceed with the measurements, report per system, and stop and report at the first failed condition.

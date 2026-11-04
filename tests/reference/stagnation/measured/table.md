# Measured sweep

All energies are Ha; residuals are Bohr^-3/2. Full windows/audit in JSON.

| Model | Policy | Arm | State | Iterations | Final residual | Certified energy | P11D triggered rules | Recovery ΔE | Sampled gap | Validated intervention |
|---|---|---|---|---:|---:|---:|---|---:|---:|---|
| cubic | reference | reference | Converged | 20 | 7.54749250509965e-11 | -1.4073901394454906 | none | None | 0.4963928002684176 | None |
| cubic | pulay-0.3-davidson | baseline | Converged | 12 | 1.8141746789098617e-09 | -1.4073901394454946 | none | None | 0.49639280029854305 | None |
| cubic | linear-1.0-davidson | baseline | Converged | 12 | 1.0322437225174038e-08 | -1.4073901394454928 | none | None | 0.49639279912058143 | None |
| cubic | linear-0.1-davidson | baseline | Converged | 119 | 8.307825684119176e-08 | -1.4073901394454875 | none | None | 0.4963927510045999 | None |
| cubic | pulay-1.0-davidson | baseline | Converged | 7 | 2.3175641641598461e-10 | -1.4073901394454942 | none | None | 0.4963928003411817 | None |
| cubic | kerker-0.3-davidson | baseline | Converged | 15 | 1.833739508968652e-08 | -1.407390139445491 | none | None | 0.4963928033989209 | None |
| cubic | linear-0.5-hybrid | baseline | Converged | 20 | 4.7659103920938736e-08 | -1.4073901394454904 | none | None | 0.4963927723667406 | None |
| cubic | kerker-0.3-hybrid | baseline | Converged | 15 | 1.833739508968652e-08 | -1.407390139445491 | none | None | 0.4963928033989209 | None |
| skewed | reference | reference | Converged | 19 | 5.796688160226345e-08 | -0.6441418148101661 | none | None | 0.024375866177302696 | None |
| skewed | pulay-0.3-davidson | baseline | Converged | 15 | 9.556911838670119e-09 | -0.644141814810181 | none | None | 0.024375907342568115 | None |
| skewed | linear-1.0-davidson | baseline | Exhausted | 200 | 0.0862306839476192 | None | Oscillation, Stagnation | None | 0.06568253953041847 | None |
| skewed | linear-1.0-davidson | intervention | Converged | 39 | 6.667116061439711e-08 | -0.6441418148101419 | Stagnation | 2.4202861936828413e-14 | 0.02437591106872128 | True |
| skewed | linear-0.1-davidson | baseline | Converged | 131 | 8.434413504314962e-08 | -0.6441418148101734 | none | None | 0.024375888477022323 | None |
| skewed | pulay-1.0-davidson | baseline | Converged | 11 | 1.446898660867456e-08 | -0.6441418148101787 | none | None | 0.02437590963216274 | None |
| skewed | kerker-0.3-davidson | baseline | Converged | 16 | 3.752657633695162e-08 | -0.6441418148101703 | none | None | 0.024375894732190595 | None |
| skewed | linear-0.5-hybrid | baseline | Converged | 22 | 5.1039973804430426e-08 | -0.644141814810177 | none | None | 0.02437589763113246 | None |
| skewed | kerker-0.3-hybrid | baseline | Converged | 16 | 3.752657633695162e-08 | -0.6441418148101703 | none | None | 0.024375894732190595 | None |
| elongated | reference | reference | Converged | 15 | 3.270324073359441e-08 | -2.873808223242774 | none | None | 0.3637878479983832 | None |
| elongated | pulay-0.3-davidson | baseline | Converged | 11 | 2.5543503861801546e-09 | -2.8738082232427775 | none | None | 0.36378786476775815 | None |
| elongated | linear-1.0-davidson | baseline | Converged | 11 | 4.456897111878379e-09 | -2.873808223242771 | none | None | 0.3637878698638542 | None |
| elongated | linear-0.1-davidson | baseline | Converged | 121 | 8.85992105123864e-08 | -2.8738082232427558 | none | None | 0.36378772117198693 | None |
| elongated | pulay-1.0-davidson | baseline | Converged | 7 | 1.613755525664115e-09 | -2.8738082232427744 | none | None | 0.3637878638955861 | None |
| elongated | kerker-0.3-davidson | baseline | Converged | 13 | 7.113599481946178e-09 | -2.873808223242773 | none | None | 0.3637878602725013 | None |
| elongated | linear-0.5-hybrid | baseline | Converged | 21 | 3.445113136224635e-08 | -2.8738082232427677 | none | None | 0.3637878024938844 | None |
| elongated | kerker-0.3-hybrid | baseline | Converged | 13 | 7.113599481946178e-09 | -2.873808223242773 | none | None | 0.3637878602725013 | None |
| slab | reference | reference | Exhausted | 1000 | 0.01331797371169278 | None | Oscillation, Stagnation | None | 0.008741195793417672 | None |
| slab | pulay-0.3-davidson | baseline | Exhausted | 200 | 0.0530841058872408 | None | Oscillation, Stagnation | None | 0.003268562955275653 | None |
| slab | pulay-0.3-davidson | intervention | Exhausted | 200 | 0.01596379424805937 | None | Oscillation, Stagnation | None | 0.025959862003517478 | False |
| slab | linear-1.0-davidson | baseline | Exhausted | 200 | 0.10778278863689246 | None | Oscillation, Stagnation | None | 0.09603129527773102 | None |
| slab | linear-1.0-davidson | intervention | Exhausted | 200 | 0.06229529457043223 | None | Oscillation, Stagnation | None | 0.004885319508185423 | False |
| slab | linear-0.1-davidson | baseline | Exhausted | 200 | 0.014179311147416943 | None | none | None | 0.01514215656801153 | None |
| slab | pulay-1.0-davidson | baseline | Exhausted | 200 | 0.05541099591357806 | None | Oscillation, Stagnation | None | 0.00726350122260562 | None |
| slab | pulay-1.0-davidson | intervention | Exhausted | 200 | 0.05420290447994257 | None | Oscillation, Stagnation | None | 0.01256701428485961 | False |
| slab | kerker-0.3-davidson | baseline | Exhausted | 200 | 0.03388231327311075 | None | Oscillation, Stagnation | None | 0.161359893917019 | None |
| slab | kerker-0.3-davidson | intervention | Exhausted | 200 | 0.05384116585816019 | None | Oscillation, Stagnation | None | 0.0143170356809343 | False |
| slab | linear-0.5-hybrid | baseline | Exhausted | 200 | 0.04845235496254963 | None | none | None | 0.28331229744868097 | None |
| slab | kerker-0.3-hybrid | baseline | Exhausted | 200 | 0.03388231327311075 | None | Oscillation, Stagnation | None | 0.161359893917019 | None |
| slab | kerker-0.3-hybrid | intervention | Exhausted | 200 | 0.05384116585816019 | None | Oscillation, Stagnation | None | 0.0143170356809343 | False |
| four-ion | reference | reference | Converged | 10 | 5.752569113677692e-08 | -3.9188912853124367 | none | None | 0.32043947318926 | None |
| four-ion | pulay-0.3-davidson | baseline | Converged | 10 | 1.407258823447631e-09 | -3.918891285312473 | none | None | 0.32043942904218914 | None |
| four-ion | linear-1.0-davidson | baseline | Converged | 7 | 3.3976634833865223e-09 | -3.9188912853124727 | none | None | 0.320439439614035 | None |
| four-ion | linear-0.1-davidson | baseline | Converged | 130 | 4.3225779741112673e-08 | -3.918891285312472 | none | None | 0.32043944057781815 | None |
| four-ion | pulay-1.0-davidson | baseline | Converged | 6 | 4.994789986290197e-09 | -3.9188912853124704 | none | None | 0.3204394106647379 | None |
| four-ion | kerker-0.3-davidson | baseline | Converged | 11 | 3.1009611667631114e-08 | -3.918891285312462 | none | None | 0.32043945880929425 | None |
| four-ion | linear-0.5-hybrid | baseline | Converged | 21 | 2.8035132412744347e-08 | -3.918891285312466 | none | None | 0.320439426279745 | None |
| four-ion | kerker-0.3-hybrid | baseline | Converged | 11 | 3.1009611667631114e-08 | -3.918891285312462 | none | None | 0.32043945880929425 | None |
| soft-mesh | reference | reference | Converged | 14 | 3.6752377987742618e-09 | -0.7090135061199141 | none | None | 0.004163613860241472 | None |
| soft-mesh | pulay-0.3-davidson | baseline | Converged | 11 | 3.0299395113274872e-09 | -0.7090135061199159 | none | None | 0.004163613396029095 | None |
| soft-mesh | linear-1.0-davidson | baseline | Converged | 8 | 1.963375398020374e-09 | -0.7090135061199156 | none | None | 0.004163612635282726 | None |
| soft-mesh | linear-0.1-davidson | baseline | Converged | 115 | 8.62717880532194e-08 | -0.7090135061199074 | none | None | 0.004163616281139038 | None |
| soft-mesh | pulay-1.0-davidson | baseline | Converged | 7 | 2.921374893108072e-09 | -0.7090135061199148 | none | None | 0.004163611793906737 | None |
| soft-mesh | kerker-0.3-davidson | baseline | Converged | 13 | 2.3713691448898028e-08 | -0.7090135061199145 | none | None | 0.0041636181043317524 | None |
| soft-mesh | linear-0.5-hybrid | baseline | Converged | 20 | 3.1409766253681285e-08 | -0.7090135061199141 | none | None | 0.004163609557676265 | None |
| soft-mesh | kerker-0.3-hybrid | baseline | Converged | 13 | 2.3713691448898028e-08 | -0.7090135061199145 | none | None | 0.0041636181043317524 | None |
| elongated-offset | reference | reference | Converged | 14 | 4.448336119512794e-08 | -2.825226251949695 | none | None | 0.2961680204705635 | None |
| elongated-offset | pulay-0.3-davidson | baseline | Converged | 13 | 5.861093972154195e-09 | -2.825226251949706 | none | None | 0.29616809224229157 | None |
| elongated-offset | linear-1.0-davidson | baseline | Converged | 40 | 1.3316574775892423e-08 | -2.8252262519496902 | none | None | 0.29616802139378584 | None |
| elongated-offset | linear-0.1-davidson | baseline | Converged | 131 | 6.23658452854277e-08 | -2.825226251949669 | none | None | 0.296167504421927 | None |
| elongated-offset | pulay-1.0-davidson | baseline | Converged | 10 | 2.6282064062115457e-10 | -2.825226251949701 | none | None | 0.29616809908047587 | None |
| elongated-offset | kerker-0.3-davidson | baseline | Converged | 24 | 7.293031856885013e-08 | -2.825226251949232 | none | None | 0.29616803665577157 | None |
| elongated-offset | linear-0.5-hybrid | baseline | Converged | 22 | 3.6521092475681576e-08 | -2.8252262519496933 | none | None | 0.2961677734099229 | None |
| elongated-offset | kerker-0.3-hybrid | baseline | Converged | 24 | 7.293031856885013e-08 | -2.825226251949232 | none | None | 0.29616803665577157 | None |
| slab-offset | reference | reference | Converged | 32 | 4.786518210819649e-08 | -0.721558424816071 | none | None | 0.03815640571276868 | None |
| slab-offset | pulay-0.3-davidson | baseline | Converged | 27 | 2.4804105159512426e-08 | -0.7215584248160546 | none | None | 0.03815640460413673 | None |
| slab-offset | linear-1.0-davidson | baseline | Exhausted | 200 | 0.11227601418830659 | None | Oscillation, Stagnation | None | 0.07989431112433508 | None |
| slab-offset | linear-1.0-davidson | intervention | Converged | 31 | 2.1375552115903586e-08 | -0.7215584248160892 | Stagnation | 1.8207657603852567e-14 | 0.03815640113315731 | True |
| slab-offset | linear-0.1-davidson | baseline | Converged | 127 | 8.329839069333561e-08 | -0.7215584248160248 | none | None | 0.03815638936749298 | None |
| slab-offset | pulay-1.0-davidson | baseline | Converged | 29 | 7.827158247853595e-09 | -0.7215584248160883 | none | None | 0.03815640092574156 | None |
| slab-offset | kerker-0.3-davidson | baseline | Converged | 45 | 1.5219251026158134e-08 | -0.7215584248160777 | Oscillation | None | 0.0381564013987913 | None |
| slab-offset | kerker-0.3-davidson | intervention | Converged | 37 | 2.4362833399216704e-08 | -0.7215584248160852 | Oscillation | 1.4210854715202004e-14 | 0.03815640432442624 | True |
| slab-offset | linear-0.5-hybrid | baseline | Exhausted | 200 | 0.042297496391354564 | None | Oscillation, Stagnation | None | 0.13249773903893206 | None |
| slab-offset | linear-0.5-hybrid | intervention | Converged | 33 | 3.1922254864187616e-08 | -0.7215584248160788 | Oscillation | 7.771561172376096e-15 | 0.03815640346930549 | True |
| slab-offset | kerker-0.3-hybrid | baseline | Converged | 45 | 1.5219251026158134e-08 | -0.7215584248160777 | Oscillation | None | 0.0381564013987913 | None |
| slab-offset | kerker-0.3-hybrid | intervention | Converged | 37 | 2.4362833399216704e-08 | -0.7215584248160852 | Oscillation | 1.4210854715202004e-14 | 0.03815640432442624 | True |
| skewed-mesh | reference | reference | Failed | 15 | 1.5588591028108474e-08 | None | none | None | None | None |
| skewed-mesh | pulay-0.3-davidson | baseline | Failed | 12 | 1.5360122300415325e-08 | None | none | None | None | None |
| skewed-mesh | linear-1.0-davidson | baseline | Failed | 71 | 7.207965912465093e-08 | None | none | None | None | None |
| skewed-mesh | linear-0.1-davidson | baseline | Failed | 125 | 8.683718036023274e-08 | None | none | None | None | None |
| skewed-mesh | pulay-1.0-davidson | baseline | Failed | 9 | 6.36832298997869e-09 | None | none | None | None | None |
| skewed-mesh | kerker-0.3-davidson | baseline | Failed | 13 | 1.4326210897423628e-08 | None | none | None | None | None |
| skewed-mesh | linear-0.5-hybrid | baseline | Failed | 21 | 4.8924587278742774e-08 | None | none | None | None | None |
| skewed-mesh | kerker-0.3-hybrid | baseline | Failed | 13 | 1.4326210897423628e-08 | None | none | None | None | None |
| four-ion-offset | reference | reference | Converged | 13 | 7.619146524401991e-09 | -3.82364672359066 | none | None | 0.3213402646640845 | None |
| four-ion-offset | pulay-0.3-davidson | baseline | Converged | 11 | 9.476564222821287e-09 | -3.823646723590664 | none | None | 0.3213402659741077 | None |
| four-ion-offset | linear-1.0-davidson | baseline | Converged | 13 | 3.12613892762473e-08 | -3.8236467235906058 | none | None | 0.3213402740851299 | None |
| four-ion-offset | linear-0.1-davidson | baseline | Converged | 127 | 4.990196417755386e-08 | -3.823646723590652 | none | None | 0.32134025832510016 | None |
| four-ion-offset | pulay-1.0-davidson | baseline | Converged | 8 | 3.218362764997769e-09 | -3.8236467235906626 | none | None | 0.3213402622513109 | None |
| four-ion-offset | kerker-0.3-davidson | baseline | Converged | 16 | 1.3159891662031193e-08 | -3.823646723590656 | none | None | 0.32134030097115474 | None |
| four-ion-offset | linear-0.5-hybrid | baseline | Converged | 20 | 4.8623698334023455e-08 | -3.8236467235906497 | none | None | 0.3213402410904027 | None |
| four-ion-offset | kerker-0.3-hybrid | baseline | Converged | 16 | 1.3159891662031193e-08 | -3.823646723590656 | none | None | 0.32134030097115474 | None |
| near-metal | reference | reference | Failed | 36 | 5.254660551095941e-08 | None | none | None | None | None |
| near-metal | pulay-0.3-davidson | baseline | Failed | 27 | 3.847083286078864e-08 | None | none | None | None | None |
| near-metal | linear-1.0-davidson | baseline | Exhausted | 200 | 0.034281493568714315 | None | Oscillation, Stagnation | None | -0.08543515220648518 | None |
| near-metal | linear-1.0-davidson | intervention | Failed | 102 | 4.0847943183810035e-08 | None | Stagnation | None | None | False |
| near-metal | linear-0.1-davidson | baseline | Failed | 155 | 9.224112305043959e-08 | None | none | None | None | None |
| near-metal | pulay-1.0-davidson | baseline | Failed | 23 | 1.3665726726416382e-09 | None | none | None | None | None |
| near-metal | kerker-0.3-davidson | baseline | Failed | 30 | 2.082164533293832e-08 | None | none | None | None | None |
| near-metal | linear-0.5-hybrid | baseline | Exhausted | 200 | 0.016757006760398926 | None | Oscillation, Stagnation | None | -0.08717156055487502 | None |
| near-metal | linear-0.5-hybrid | intervention | Failed | 70 | 5.1136094988817606e-08 | None | Stagnation | None | None | False |
| near-metal | kerker-0.3-hybrid | baseline | Failed | 30 | 2.082164533293832e-08 | None | none | None | None | None |

# Performance benchmark

Development workstation:

* Julia 1.12.6
* 28 Julia threads
* initial H-ENA altitude: 600 km
* initial velocity: [-400, 0, 0] km/s
* solar minimum, Ls = 0 degrees
* MGITM cold atmosphere plus MAMPS hot O
* maximum transport step: 1 km

| Particles | Threads | Elapsed (s) | Throughput (particles/s) |
| ---: | ---: | ---: | ---: |
| 1,000,000 | 28 | 131.042 | 7,631 |

The 1,000,000-particle run produced:

* mean final energy: 177.210 eV
* mean collision count: 1094.741
* mean transport step count: 4253.976
* stopped below minimum energy: 420,973
* escaped above the maximum altitude: 423,793
* reached the maximum collision count: 155,234

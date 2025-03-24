# Iceplume Tracer Stability

In running this model, we encountered numerical instabilities associated with the treatment of tracer quantities on days that transitioned from low to exceptionally high discharge. This instability was not observed in the temperature and salinity fields - only the passive tracers - and it was not alleviated by reducing timestepping or changing the advective scheme. To circumvent this issue, we applied bounds to the tracer fields in the plume calculation. We document this here:

| Tracer Number | Tracer Quantity in L2_Disko_Bay configuration | Bounds |
|---------------|-----------------------------------------------|--------|
| 1             | DIC                                           | -20 to 5000 uM |
| 2             | NO3                                           | 0 to 30 uM |
| 3             | NO2                                           | 0 to 30 uM |
| 4             | NH4                                           | 0 to 30 uM |
| 5             | PO4                                           | 0 to 30 uM |
| 6             | FeT                                           | 0 to 3 uM |
| 7             | SiO2                                           | 0 to 30 uM |
| 8             | DOC                                           | 0 to 300 uM |
| 9             | DON                                           | 0 to 30 uM |
| 10            | DOP                                           | 0 to 30 uM |
| 11            | DOFe                                           | 0 to 0.3 uM |
| 12            | POC                                           | 0 to 30 uM |
| 13            | PON                                           | 0 to 30 uM |
| 14            | POP                                           | 0 to 30 uM |
| 15            | POFe                                           | 0 to 0.3 uM |
| 16            | POSi                                           | 0 to 30 uM |
| 17            | PIC                                           | 0 to 30 uM |
| 18            | Alk                                           | 0 to 5000 meq/m3 |
| 19            | O2                                           | 0 to 1000 uM |
| 20            | c01                                           | -1 to 50 uM |
| 21            | c02                                           | -1 to 50 uM |
| 22            | c03                                           | -1 to 50 uM |
| 23            | c04                                           | -1 to 50 uM |
| 24            | c05                                           | -1 to 50 uM |
| 25            | c06                                           | -1 to 50 uM |
| 26            | c07                                           | -1 to 50 uM |
| 27            | Chl01                                           | 0 to 30 mg/m3 |
| 28            | Chl02                                           | 0 to 30 mg/m3 |
| 29            | Chl03                                           | 0 to 30 mg/m3 |
| 30            | Chl04                                           | 0 to 30 mg/m3 |
| 31            | Chl05                                           | 0 to 30 mg/m3 |




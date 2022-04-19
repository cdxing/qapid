# qapid
QA/PID for star pico dst analysis, this package is part of a STAR Flow analysis

## Prerequisite
StRoot/StPicoEvent 

## Steps
starver XXXX # the corresponding version of the dataset, here we use the year2018 3GeV dataset as an example

```
starver SL19e
cons
./run.sh

```

## Clear and create new folders for scheduler jobs
```
./clear.sh
```

## Submit scheduler jobs
```
star-submit submit.xml
```


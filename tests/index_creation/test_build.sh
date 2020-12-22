#!/bin/bash

../../HowDeSBT_strains/howdesbt build --howde --tree=union.sbt --outtree=howde.sbt

diff ./howde_test.sbt ./howde.sbt

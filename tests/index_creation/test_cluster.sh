#!/bin/bash

../../HowDeSBT_strains/howdesbt cluster --list=leafname_merge --tree=union.sbt --nodename=node{number} --cull

diff ./union_test.sbt ./union.sbt

#!/bin/bash

echo "Build index ..."
howdesbt cluster --list=leafname --tree=union.sbt --nodename=node{number} --cull
howdesbt build --howde --tree=union.sbt --outtree=howde.sbt
ls | grep -Pv 'detbrief.rrr.' | grep '.bf' | xargs rm --
echo "Done"

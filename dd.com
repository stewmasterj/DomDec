        1540 # max ourNpts across all images
   20.0000000       20.0000000       20.0000000     # global BoxSize
           1        1240        1540           3 # imID, myNpts, ourNpts, NNcells
    2    3    4 # NeighCellIDs
 bulkminmax:           1         972
           3 # NumPermutationsSendCells
sendmin:       973      1081      1225
sendmax:      1080      1224      1240
           5 # NumPermutationsGetCells
getmin:      1241      1349      1365      1509      1525
getmax:      1348      1364      1508      1524      1540
sourceCellIDs:    2    2    3    3    4
sourcemin:       973      1081       973      1117       973
sourcemax:      1080      1096      1116      1132       988
           2        1240        1540           3 # imID, myNpts, ourNpts, NNcells
    1    3    4 # NeighCellIDs
 bulkminmax:           1         972
           3 # NumPermutationsSendCells
sendmin:       973      1081      1097
sendmax:      1080      1096      1240
           5 # NumPermutationsGetCells
getmin:      1241      1349      1365      1381      1397
getmax:      1348      1364      1380      1396      1540
sourceCellIDs:    1    1    3    4    4
sourcemin:       973      1225      1117       973       989
sourcemax:      1080      1240      1132       988      1132
           3        1240        1540           3 # imID, myNpts, ourNpts, NNcells
    1    2    4 # NeighCellIDs
 bulkminmax:           1         972
           3 # NumPermutationsSendCells
sendmin:       973      1117      1133
sendmax:      1116      1132      1240
           5 # NumPermutationsGetCells
getmin:      1241      1385      1401      1417      1433
getmax:      1384      1400      1416      1432      1540
sourceCellIDs:    1    1    2    4    4
sourcemin:      1081      1225      1081       973      1133
sourcemax:      1224      1240      1096       988      1240
           4        1240        1540           3 # imID, myNpts, ourNpts, NNcells
    1    2    3 # NeighCellIDs
 bulkminmax:           1         972
           3 # NumPermutationsSendCells
sendmin:       973       989      1133
sendmax:       988      1132      1240
           5 # NumPermutationsGetCells
getmin:      1241      1257      1273      1417      1433
getmax:      1256      1272      1416      1432      1540
sourceCellIDs:    1    2    2    3    3
sourcemin:      1225      1081      1097      1117      1133
sourcemax:      1240      1096      1240      1132      1240

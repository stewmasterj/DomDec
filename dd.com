        1650 # max ourNpts across all images
   10.0000000       10.0000000       1.00000000     # global BoxSize
           1         766        1105           4 # imID, myNpts, ourNpts, NNcells
    2    3    6    7 # NeighCellIDs
 bulkminmax:           1         480
           5 # NumPermutationsSendCells
sendmin:       481       635       658       719       731
sendmax:       634       657       718       730       766
           8 # NumPermutationsGetCells
getmin:       767       794       829       975      1005      1034      1047      1071
getmax:       793       828       974      1004      1033      1046      1070      1105
sourceCellIDs:    2    2    3    3    6    6    6    7
sourcemin:       316       343       361       635       126       155       168       186
sourcemax:       342       377       506       664       154       167       191       220
           2         766        1190           2 # imID, myNpts, ourNpts, NNcells
    1    6 # NeighCellIDs
 bulkminmax:           1         315
           3 # NumPermutationsSendCells
sendmin:       316       343       378
sendmax:       342       377       766
           4 # NumPermutationsGetCells
getmin:       767       790       851       880
getmax:       789       850       879      1190
sourceCellIDs:    1    1    6    6
sourcemin:       635       658       126       192
sourcemax:       657       718       154       502
           3         766        1360           5 # imID, myNpts, ourNpts, NNcells
    1    4    6    7    8 # NeighCellIDs
 bulkminmax:           1         360
           7 # NumPermutationsSendCells
sendmin:       361       507       635       665       692       739       750
sendmax:       506       634       664       691       738       749       766
          11 # NumPermutationsGetCells
getmin:       767       921       957      1129      1163      1187      1222      1245      1278      1298      1312
getmax:       920       956      1128      1162      1186      1221      1244      1277      1297      1311      1360
sourceCellIDs:    1    1    4    4    6    7    7    7    8    8    8
sourcemin:       481       731       301       633       168       186       221       244       391       411       425
sourcemax:       634       766       472       666       191       220       243       276       410       424       473
           4         766        1270           3 # imID, myNpts, ourNpts, NNcells
    3    5    8 # NeighCellIDs
 bulkminmax:           1         300
           5 # NumPermutationsSendCells
sendmin:       301       473       633       667       736
sendmax:       472       632       666       735       766
           9 # NumPermutationsGetCells
getmin:       767       895       906       923      1063      1090      1104      1153      1209
getmax:       894       905       922      1062      1089      1103      1152      1208      1270
sourceCellIDs:    3    3    3    5    5    8    8    8    8
sourcemin:       507       739       750       571       711       411       425       474       530
sourcemax:       634       749       766       710       737       424       473       529       591
           5         765        1020           2 # imID, myNpts, ourNpts, NNcells
    4    8 # NeighCellIDs
 bulkminmax:           1         570
           3 # NumPermutationsSendCells
sendmin:       571       711       738
sendmax:       710       737       765
           4 # NumPermutationsGetCells
getmin:       766       926       957      1019
getmax:       925       956      1018      1020
sourceCellIDs:    4    4    8    8
sourcemin:       473       736       530       592
sourcemax:       632       766       591       593
           6         766        1650           4 # imID, myNpts, ourNpts, NNcells
    1    2    3    7 # NeighCellIDs
 bulkminmax:           1         125
           5 # NumPermutationsSendCells
sendmin:       126       155       168       192       503
sendmax:       154       167       191       502       766
           8 # NumPermutationsGetCells
getmin:       767       828       840       876       911      1300      1330      1365
getmax:       827       839       875       910      1299      1329      1364      1650
sourceCellIDs:    1    1    1    2    2    3    7    7
sourcemin:       658       719       731       343       378       635       186       277
sourcemax:       718       730       766       377       766       664       220       562
           7         765        1410           4 # imID, myNpts, ourNpts, NNcells
    1    3    6    8 # NeighCellIDs
 bulkminmax:           1         185
           5 # NumPermutationsSendCells
sendmin:       186       221       244       277       563
sendmax:       220       243       276       562       765
          10 # NumPermutationsGetCells
getmin:       766       802       832       859       906       917       941      1205      1225      1239
getmax:       801       831       858       905       916       940      1204      1224      1238      1410
sourceCellIDs:    1    3    3    3    3    6    6    8    8    8
sourcemin:       731       635       665       692       739       168       503       391       411       594
sourcemax:       766       664       691       738       749       191       766       410       424       765
           8         765        1265           4 # imID, myNpts, ourNpts, NNcells
    3    4    5    7 # NeighCellIDs
 bulkminmax:           1         390
           7 # NumPermutationsSendCells
sendmin:       391       411       425       474       530       592       594
sendmax:       410       424       473       529       591       593       765
          10 # NumPermutationsGetCells
getmin:       766       813       824       841       875       944       975      1002      1030      1063
getmax:       812       823       840       874       943       974      1001      1029      1062      1265
sourceCellIDs:    3    3    3    4    4    4    5    5    7    7
sourcemin:       692       739       750       633       667       736       711       738       244       563
sourcemax:       738       749       766       666       735       766       737       765       276       765

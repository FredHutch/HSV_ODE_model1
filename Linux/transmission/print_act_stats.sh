#!/bin/sh
echo "No drug results..."
grep acts no_drug/all_trans.out | awk '{r++;if ($1 == "No") sum += $4; else {n++; sum+= $3;}} END {printf("%d runs, %d acts, %d transmissions, %lf avg\n",r,sum,n,sum/n);}'
echo 
echo "Regime1 results..."
grep acts regime1/all_trans.out | awk '{r++;if ($1 == "No") sum += $4; else {n++; sum+= $3;}} END {printf("%d runs, %d acts, %d transmissions, %lf avg\n",r,sum,n,sum/n);}'
echo 
echo "Regime2 results..."
grep acts regime2/all_trans.out | awk '{r++;if ($1 == "No") sum += $4; else {n++; sum+= $3;}} END {printf("%d runs, %d acts, %d transmissions, %lf avg\n",r,sum,n,sum/n);}'
echo 
echo "Regime3 results..."
grep acts regime3/all_trans.out | awk '{r++;if ($1 == "No") sum += $4; else {n++; sum+= $3;}} END {printf("%d runs, %d acts, %d transmissions, %lf avg\n",r,sum,n,sum/n);}'
echo 
echo "Regime4 results..."
grep acts regime4/all_trans.out | awk '{r++;if ($1 == "No") sum += $4; else {n++; sum+= $3;}} END {printf("%d runs, %d acts, %d transmissions, %lf avg\n",r,sum,n,sum/n);}'
echo 
echo "Regime5 results..."
grep acts regime5/all_trans.out | awk '{r++;if ($1 == "No") sum += $4; else {n++; sum+= $3;}} END {printf("%d runs, %d acts, %d transmissions, %lf avg\n",r,sum,n,sum/n);}'
echo 
echo "Regime6 results..."
grep acts regime6/all_trans.out | awk '{r++;if ($1 == "No") sum += $4; else {n++; sum+= $3;}} END {printf("%d runs, %d acts, %d transmissions, %lf avg\n",r,sum,n,sum/n);}'

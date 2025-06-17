# Running all instances of the sets E and M
./VRPSolver data/E/E-n51-k5.vrp -m 5 -M 5 -u 522 &> out/E-n51-k5.out
./VRPSolver data/E/E-n76-k7.vrp -m 7 -M 7 -u 683 &> out/E-n76-k7.out
./VRPSolver data/E/E-n76-k8.vrp -m 8 -M 8 -u 736 &> out/E-n76-k8.out
./VRPSolver data/E/E-n76-k10.vrp -m 10 -M 10 -u 831 &> out/E-n76-k10.out
./VRPSolver data/E/E-n76-k14.vrp -m 14 -M 14 -u 1022 &> out/E-n76-k14.out
./VRPSolver data/E/E-n101-k8.vrp -m 8 -M 8 -u 816 &> out/E-n101-k8.out
./VRPSolver data/E/E-n101-k14.vrp -m 14 -M 14 -u 1072 &> out/E-n101-k14.out
./VRPSolver data/M/M-n101-k10.vrp -m 10 -M 10 -u 821 &> data/M/M-n101-k10.out
./VRPSolver data/M/M-n121-k7.vrp -m 7 -M 7 -u 1035 &> data/M/M-n121-k7.out
./VRPSolver data/M/M-n151-k12.vrp -m 12 -M 12 -u 1016 &> data/M/M-n151-k12.out
# M-n200-k16 takes several hours     
#./VRPSolver data/M/M-n200-k16.vrp -m 16 -M 16 -u 1279 &> data/M/M-n200-k16.out
./VRPSolver data/M/M-n200-k17.vrp -m 17 -M 17 -u 1276 &> data/M/M-n200-k17.out

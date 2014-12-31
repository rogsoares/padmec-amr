#simple script for the sequence :
#make ,go to simulator's dir,make,(conditional)./simulator.exe
clear
clear
clear

echo"cleanning->multigridlib dir..."
make cleanup
echo "make->multigrid"
echo ""
echo ""
make
echo ""

#simulator directory
cd ../../../prototypes/SIMULATORPROJ1
echo "cleanning->simulator dir..."
make cleanup
echo "make->simulator"
echo ""
echo ""

make

echo  "run the simulator now? "
read value

	if [ "$value" -eq "1" ];
        then
	    echo "running simulator.exe..."
	    sleep 1
		./simulator.exe
	fi


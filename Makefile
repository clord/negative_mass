sim: sim.cc 
	g++ -x objective-c++ -framework Foundation -framework AppKit -framework QTKit sim.cc -o sim

clean:
	rm sim

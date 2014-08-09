include ../../../environment

all: MainSeq  MainLan 

clean:
	rm -f MainLan MainSeq MainWan *.o *% *~


MainSeq: PSO.req.o PSO.pro.o MainSeq.o
	$(CXX) $(LDFLAGS) $^ $(LOADLIBES) $(CPPFLAGS) -o $@

MainLan: PSO.req.o PSO.pro.o MainLan.o
	$(CXX) $(LDFLAGS) $^ $(LOADLIBES) $(CPPFLAGS) -o $@

#MainWan: PSO.req.o PSO.pro.o MainWan.o
#	$(CXX) $(LDFLAGS) $^ $(LOADLIBES) $(CPPFLAGS) -o $@


LAN:
#	$(RUN)  -v  -p4pg pgfileLan MainLan
#	$(RUN)  -v  -np 4 MainLan PSO.cfg ../../../ProblemInstances/ONEMAX-instances/onemax10.txt res/om10.newga.seq.txt
	$(RUN)  -v  -np 2 MainLan PSO.cfg ../../../ProblemInstances/ONEMAX-instances/onemax10.txt res/om10.newga.seq.txt
#WAN:
#	$(RUN)  -v  -p4pg pgfileWan MainWan
SEQ:
	./MainSeq PSO.cfg ../../../ProblemInstances/ONEMAX-instances/onemax10.txt res/om10.newga.seq.txt

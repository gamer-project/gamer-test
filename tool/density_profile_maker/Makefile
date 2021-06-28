

CC:=g++
exe:=main
obj:=main.o

INCDIR += -I$(GSL_DIR)/include 
LIBS += -L$(GSL_DIR)/lib 

all:$(obj)
	$(CC) -o $(exe) $(obj)  $(LIBS) 
%.o:%.c $(INCS)
	$(CC) -I$(INCDIR) -c $^ -o $@ 

.PHONY:clean
clean:
	rm -rf $(obj) $(exe)


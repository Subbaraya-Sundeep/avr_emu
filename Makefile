TARGET := emulator
C_SRCS := $(wildcard *.c)
C_OBJS := ${C_SRCS:.c=.o}

INCLUDE_DIRS := ./include
CC := gcc

CPPFLAGS += $(foreach includedir,$(INCLUDE_DIRS),-I$(includedir))

all: $(TARGET)

$(TARGET): $(C_OBJS)
	$(CC) $(C_OBJS) -o $(TARGET) -lpthread

clean :
	rm -f *.o
	rm -f emulator

# Compiler and Flags
CC = g++
CFLAGS = -m64 -O3 -shared -fPIC -pg
BOOST_PATH = /home/amk7r/colocation_mining/arpan_codes/Codes/boost/

# Source and Target
SRCS = functions.cpp
TARGET = functions.so

# Build Rule
$(TARGET): $(SRCS)
	$(CC) $(CFLAGS) -I$(BOOST_PATH) -o $@ $(SRCS)

# Clean Rule
clean:
	rm -f $(TARGET)

# Run Rule
run: $(TARGET)
	./$(TARGET)

CC = clang++
CPPFLAGS = 
OUTPUT = page_rank
SRC_DIR = src
SRC = $(wildcard $(SRC_DIR)/*.cpp)
OBJ = $(SRC:$(SRC_DIR)/%.cpp=$(SRC_DIR)/%.o)


$(SRC_DIR)/%.o : $(SRC_DIR)/%.cpp
	$(CC)  -o $@ -c  $< $(CPPFLAGS)  

$(OUTPUT) : $(OBJ)
	$(CC)  -o $@ $^ $(CPPFLAGS)


clean :
	rm -f $(OBJ) $(OUTPUT)
CC := gcc
CFLAGS := -Wall -O3 -fopenmp
OBJFLAGS := $(CFLAGS) -c 
INC := -I include

SRCDIR := src
BUILDDIR := build
TARGET := bin/theta-method

SOURCES := $(shell find $(SRCDIR) -type f -name '*.c')
OBJECTS := $(patsubst $(SRCDIR)/%, $(BUILDDIR)/%, $(SOURCES:.c=.o))
DEPS := $(OBJECTS:%.o=%.d)

.SECONDEXPANSION:
.PHONY: clean

all: $(TARGET)

$(TARGET): $(OBJECTS)
	@echo -e "Linking"
	@mkdir -p $$(dirname $@)
	$(CC) $(INC) -fopenmp $^ -o $@

$(BUILDDIR)/%.o: $(SRCDIR)/%.c
	@mkdir -p $$(dirname $@)
	@echo -e "Building $@"
	$(CC) $(INC) $(OBJFLAGS) -o $@ -MMD -MP -MF ${@:.o=.d} $<

clean:
	rm -rf $(BUILDDIR) $(TARGET) vgcore* core log
	rmdir bin

-include $(DEPS)


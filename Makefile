
# Makefile for Bulletproof CUDA implementation
NVCC = nvcc
NVCC_FLAGS = -arch=sm_80 -O3 -diag-suppress 177 -Xcompiler "-Wno-deprecated-declarations" -lcrypto -lssl

# Source files - organized by module
CORE_SOURCES = curve25519_ops.cu
VECTOR_SOURCES = bulletproof_vectors.cu
CHALLENGE_SOURCES = bulletproof_challenge.cu
PROTOCOL_SOURCES = bulletproof_range_proof.cu
TEST_SOURCES = complete_bulletproof_test.cu

# CUDA optimization sources
CUDA_OPT_SOURCES = cuda_bulletproof_kernels.cu cuda_inner_product.cu cuda_field_ops.cu cuda_range_proof_verify.cu

# Combine all sources
SOURCES = $(CORE_SOURCES) $(VECTOR_SOURCES) $(CHALLENGE_SOURCES) $(PROTOCOL_SOURCES) $(TEST_SOURCES)
CUDA_SOURCES = $(CUDA_OPT_SOURCES)

# Target file names for .o files
OBJECTS = $(SOURCES:.cu=.o)
CUDA_OBJECTS = $(CUDA_OPT_SOURCES:.cu=.o)

# Output binary
TARGET = cuda_bulletproof_test

# Build rule
all: $(TARGET)

# Compile each .cu file separately
%.o: %.cu
	$(NVCC) -c $(NVCC_FLAGS) $< -o $@

# Special rule for CUDA files to include device_curve25519_ops.cuh
$(CUDA_OBJECTS): %.o: %.cu device_curve25519_ops.cuh
	$(NVCC) -c $(NVCC_FLAGS) $< -o $@

# Link the object files
$(TARGET): $(OBJECTS) $(CUDA_OBJECTS)
	$(NVCC) $(NVCC_FLAGS) $(OBJECTS) $(CUDA_OBJECTS) -o $(TARGET)

# Clean rule
clean:
	rm -f $(TARGET) $(OBJECTS) $(CUDA_OBJECTS)

run: $(TARGET)
	./$(TARGET)

# Debugging targets
debug_flags = -g -G -O0
debug: NVCC_FLAGS += $(debug_flags)
debug: all

# Profiling target
profile: NVCC_FLAGS += -lineinfo
profile: all
	nvprof --metrics all ./$(TARGET)

# Test with different optimization levels
test_O0: NVCC_FLAGS += -O0
test_O0: all
	./$(TARGET)

test_O2: NVCC_FLAGS += -O2
test_O2: all
	./$(TARGET)

# Run with memory checker
memcheck: debug
	cuda-memcheck ./$(TARGET)

# Performance benchmarking
benchmark: all
	@echo "Running performance benchmarks..."
	./$(TARGET) --benchmark

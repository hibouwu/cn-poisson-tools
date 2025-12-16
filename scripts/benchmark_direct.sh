#!/bin/bash

# Compile the project
echo "Compiling..."
cd "$(dirname "$0")/.." || exit
make

# Output file
OUTPUT_FILE="benchmark_results.txt"
echo "Running benchmarks... Results will be saved to $OUTPUT_FILE"
echo "Method,Size,Time(ms)" > "$OUTPUT_FILE"

# Define sizes to test
SIZES=(100 1000 10000 100000)

# Define methods: 0=TRF (LAPACK), 1=TRI (Custom), 2=SV (LAPACK Driver)
METHODS=(0 1 2)

for size in "${SIZES[@]}"; do
    for method in "${METHODS[@]}"; do
        echo "Running Method $method with N=$size..."
        
        # Run the program and capture the output
        # We start looking for the line "Execution time..."
        result=$(./bin/tpPoisson1D_direct "$method" "$size")
        
        # Extract time using grep and awk (assuming format: "Execution time (IMPLEM=X, N=Y): Z ms")
        time_ms=$(echo "$result" | grep "Execution time" | awk '{print $(NF-1)}')
        
        if [ -z "$time_ms" ]; then
            time_ms="Error"
        fi
        
        echo "$method,$size,$time_ms" >> "$OUTPUT_FILE"
    done
done

echo "Benchmark complete."
cat "$OUTPUT_FILE"

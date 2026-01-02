#!/bin/bash

# Compile the project
echo "Compiling..."
cd "$(dirname "$0")/.." || exit
make

# Output file
OUTPUT_FILE="benchmark_results_iter.txt"
echo "Running iterative benchmarks... Results will be saved to $OUTPUT_FILE"
echo "Method,Size,Time(ms),RelRes,Iterations" > "$OUTPUT_FILE"

# Define sizes to test
SIZES=(10 100 1000)

# Define methods: 0=ALPHA (Richardson), 1=JAC (Jacobi), 2=GS (Gauss-Seidel), 3=CSR, 4=CSC
METHODS=(0 1 2 3 4)

for size in "${SIZES[@]}"; do
    for method in "${METHODS[@]}"; do
        echo "Running Iterative Method $method with N=$size..."
        
        result=$(./bin/tpPoisson1D_iter "$method" "$size")
        
        # Extract metrics
        time_ms=$(echo "$result" | grep "Execution time" | awk '{print $(NF-1)}')
        nb_ite=$(echo "$result" | grep "Nb iterations" | awk '{print $NF}')
        rel_res=$(echo "$result" | grep "relres =" | awk '{print $NF}')
        
        if [ -z "$time_ms" ]; then time_ms="Error"; fi
        if [ -z "$nb_ite" ]; then nb_ite="Error"; fi
        if [ -z "$rel_res" ]; then rel_res="Error"; fi

        echo "$method,$size,$time_ms,$rel_res,$nb_ite" >> "$OUTPUT_FILE"
    done
done

echo "Benchmark complete."
cat "$OUTPUT_FILE"

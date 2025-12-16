import matplotlib.pyplot as plt
import pandas as pd
import sys
import os

def plot_benchmark(file_path):
    # Check if file exists
    if not os.path.exists(file_path):
        print(f"Error: File '{file_path}' not found.")
        print("Please run the benchmark script first: ./scripts/benchmark_direct.sh")
        return

    # Read data
    try:
        df = pd.read_csv(file_path)
    except Exception as e:
        print(f"Error reading CSV: {e}")
        return

    # Method Labels
    method_labels = {
        0: 'LAPACK dgbtrf (Band)',
        1: 'Custom Tridiagonal',
        2: 'LAPACK dgbsv (Simple Driver)'
    }

    plt.figure(figsize=(10, 6))

    # Plot each method
    for method_id in sorted(df['Method'].unique()):
        subset = df[df['Method'] == method_id].sort_values(by='Size')
        label = method_labels.get(method_id, f'Method {method_id}')
        plt.plot(subset['Size'], subset['Time(ms)'], marker='o', label=label)

    # Configuration
    plt.title('Poisson 1D Direct Solver Performance')
    plt.xlabel('Matrix Size (N)')
    plt.ylabel('Execution Time (ms)')
    plt.grid(True, which="both", ls="-", alpha=0.5)
    plt.legend()
    
    # Log scales often allow better visualization of scaling (O(N))
    plt.xscale('log')
    plt.yscale('log')

    # Save plot
    output_file = 'benchmark_plot.png'
    plt.savefig(output_file)
    print(f"Plot saved to {output_file}")
    # plt.show()

if __name__ == "__main__":
    # Default to benchmark_results.txt in current dir, or take arg
    if len(sys.argv) > 1:
        input_file = sys.argv[1]
    else:
        input_file = 'benchmark_results.txt'
    
    plot_benchmark(input_file)

import matplotlib.pyplot as plt
import pandas as pd
import sys
import os

def plot_combined_convergence(file_label_pairs):
    plt.figure(figsize=(10, 6))
    
    colors = ['b', 'r', 'g', 'c', 'm', 'y', 'k']
    styles = ['-', '--', '-.', ':']
    
    for i, (file_path, label) in enumerate(file_label_pairs):
        if not os.path.exists(file_path):
            print(f"Warning: File '{file_path}' not found. Skipping.")
            continue

        try:
            data = pd.read_csv(file_path, header=None, names=['Residual'])
            data['Iteration'] = data.index + 1
            
            # Use distinct color/style
            color = colors[i % len(colors)]
            style = styles[i % len(styles)]
            
            plt.plot(data['Iteration'], data['Residual'], label=label, 
                     color=color, linestyle=style, linewidth=2)
            
            # Annotate the last point
            last_iter = data['Iteration'].iloc[-1]
            last_val = data['Residual'].iloc[-1]
            plt.annotate(f'{last_val:.2e}', xy=(last_iter, last_val), 
                         xytext=(10, 0), textcoords='offset points',
                         color=color, fontsize=9, va='center')
                     
        except Exception as e:
            print(f"Error reading {file_path}: {e}")
            continue

    plt.yscale('log')
    plt.title('Convergence History Comparison')
    plt.xlabel('Iteration')
    plt.ylabel('Relative Residual ||r||/||b||')
    plt.grid(True, which="both", ls="-", alpha=0.5)
    plt.legend()
    
    output_png = 'convergence_comparison.png'
    plt.savefig(output_png)
    print(f"Combined plot saved to {output_png}")

if __name__ == "__main__":
    # Expects pairs: python plot_convergence.py file1 label1 file2 label2 ...
    args = sys.argv[1:]
    if len(args) >= 2 and len(args) % 2 == 0:
        pairs = []
        for i in range(0, len(args), 2):
            pairs.append((args[i], args[i+1]))
        plot_combined_convergence(pairs)
    else:
        print("Usage: python plot_convergence.py file1 Label1 [file2 Label2 ...]")

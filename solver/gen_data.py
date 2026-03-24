"""Generate histogram data as JSON for the web visualization."""
import json
import sys
from collections import Counter
from .irred_experiment import run_one_n

def main():
    results = {}
    for n in range(3, 11):
        total_counts, irred_counts, num_valid, num_capped = run_one_n(
            n, num_seeds=200, max_solutions=5000
        )
        
        total_freq = dict(Counter(total_counts))
        irred_freq = dict(Counter(irred_counts))
        
        results[str(n)] = {
            "total": {str(k): v for k, v in sorted(total_freq.items())},
            "irreducible": {str(k): v for k, v in sorted(irred_freq.items())},
            "num_valid": num_valid,
        }
        
        print(f"n={n}: {num_valid} valid, total_keys={sorted(total_freq.keys())}, "
              f"irred_keys={sorted(irred_freq.keys())}", flush=True)
    
    with open("solver/histogram_data.json", "w") as f:
        json.dump(results, f, indent=2)
    
    print("\nData written to solver/histogram_data.json", flush=True)
    # Also print it so we can capture
    print("JSON_START")
    print(json.dumps(results))
    print("JSON_END")

if __name__ == "__main__":
    main()

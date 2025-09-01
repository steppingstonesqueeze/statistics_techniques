# alias_sampler.py
import argparse
import numpy as np
import time

def create_alias_table(probs):
    n = len(probs)
    scaled_probs = np.array(probs) * n
    alias = [0] * n
    prob = [0.0] * n
    small, large = [], []

    for i, sp in enumerate(scaled_probs):
        if sp < 1.0:
            small.append(i)
        else:
            large.append(i)

    while small and large:
        s = small.pop()
        l = large.pop()

        prob[s] = scaled_probs[s]
        alias[s] = l

        scaled_probs[l] = scaled_probs[l] - (1.0 - prob[s])
        if scaled_probs[l] < 1.0:
            small.append(l)
        else:
            large.append(l)

    for remaining in small + large:
        prob[remaining] = 1.0
        alias[remaining] = remaining

    return alias, prob

def alias_sample(alias, prob, n):
    K = len(alias)
    samples = []
    for _ in range(n):
        i = np.random.randint(K)
        if np.random.rand() < prob[i]:
            samples.append(i)
        else:
            samples.append(alias[i])
    return samples

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--k", type=int, default=10000, help="Number of categories")
    parser.add_argument("--n", type=int, default=100000, help="Number of samples to draw")
    args = parser.parse_args()

    k = args.k
    n = args.n
    probs = np.random.dirichlet(np.ones(k))

    start = time.time()
    alias, prob = create_alias_table(probs)
    build_time = time.time() - start

    start = time.time()
    samples = alias_sample(alias, prob, n)
    sample_time = time.time() - start

    print(f"Alias table built in {build_time:.3f} seconds")
    print(f"Sampled {n} values in {sample_time:.3f} seconds")
    print("First 10 samples:", samples[:10])

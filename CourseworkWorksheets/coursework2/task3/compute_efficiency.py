import sys
    
def efficiency(serialTime, parallelTime, numProcesses):
    return serialTime / (numProcesses * parallelTime)

if __name__ == "__main__":
    # Parse command line arguments
    if len(sys.argv) < 3:
        sys.exit("Usage: python {0} <inputRuntimeFile> <outputEfficiencyFile>".format(sys.argv[0]))
    inputRuntimeFile = sys.argv[1]
    outputEfficiencyFile = sys.argv[2]
    
    # Load input file
    runtimes = []
    with open(inputRuntimeFile, "r") as f:
        for line in f.readlines():
            fields = line.split()
            runtimes.append( (int(fields[0]), float(fields[1])) )
    # Find serial execution time
    serialTime = None
    for rt in runtimes:
        if rt[0] == 1:
            serialTime = rt[1]
    if not serialTime:
        sys.exit("Serial execution time not found in data file")
    # Compute efficiencies
    efficiencies = []
    for rt in runtimes:
        efficiencies.append( (rt[0], efficiency(serialTime, rt[1], rt[0])) )
    # Write speedups to a file
    with open(outputEfficiencyFile, "w") as f:
        for ef in efficiencies[:-1]: # to prevent extra empty line being written
            f.write("{0} {1}\n".format(*ef))
        f.write("{0} {1}".format(efficiencies[-1][0], efficiencies[-1][1]))

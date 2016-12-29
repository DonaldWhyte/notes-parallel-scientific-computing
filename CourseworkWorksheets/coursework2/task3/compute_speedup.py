import sys

def speedup(serialTime, parallelTime):
    return serialTime / parallelTime

if __name__ == "__main__":
    # Parse command line arguments
    if len(sys.argv) < 3:
        sys.exit("Usage: python {0} <inputRuntimeFile> <outputSpeedupFile>".format(sys.argv[0]))
    inputRuntimeFile = sys.argv[1]
    outputSpeedupFile = sys.argv[2]
    
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
    # Compute speed ups
    speedups = []
    for rt in runtimes:
        speedups.append( (rt[0], speedup(serialTime, rt[1])) )
    # Write speedups to a file
    with open(outputSpeedupFile, "w") as f:
        for su in speedups[:-1]: # to prevent extra empty line being written
            f.write("{0} {1}\n".format(*su))
        f.write("{0} {1}".format(speedups[-1][0], speedups[-1][1]))

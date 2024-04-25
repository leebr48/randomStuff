# This script calculates the mirror ratio of wout_Lee_1, wout_Lee_2, and wout_Lee_3. The min and max |B| were obtained from VMECplot

mins = [3.3142, 4.0283, 4.1522]
maxes = [8.1282, 7.1282, 7.2144]

for i in range(len(mins)):
    minUse = mins[i]
    maxUse = maxes[i]
    mirrorRat = (maxUse - minUse)/(maxUse + minUse)
    print(mirrorRat * 100)

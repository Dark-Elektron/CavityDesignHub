import matplotlib.pyplot as plt

def p3x1(v):
    pattern.append(int(v))
    if v == 1:
        return
        # print('End: ', pattern)
    else:
        if v%2 == 0:
            return p3x1(v/2)
        else:
            return p3x1(3*v + 1)

for i in range(3, 150):
    pattern = []
    p3x1(i)
    plt.plot(pattern, marker='o')

plt.show()

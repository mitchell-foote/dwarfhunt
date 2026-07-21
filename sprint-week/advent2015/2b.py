
def ribbon_ft(l, w, h):
    measurements = [l,w,h]
    min_measure =  min(measurements)
    measurements.remove(min_measure)
    second_min = min(measurements)

    return l*w*h + min_measure*2 + second_min*2

print(ribbon_ft(2,3,4))
print(ribbon_ft(1,1,10))
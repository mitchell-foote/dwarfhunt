def santa_map(path):
    x = 0
    y = 0
    arivedMap = set({(0,0)})
    for char in path:
        if char == '>':
            x += 1
        elif char == '^':
            y += 1
        elif char == '<':
            x -= 1
        else:
            y -= 1
        arivedMap.add((x,y))
    return arivedMap.__len__()

def robo_santa_map(path):
    robo_path = [x for i, x in enumerate(path) if i % 2 == 0 ]
    santa_path = [x for i, x in enumerate(path) if i % 2 != 0 ]
    return (santa_map(robo_path) - 1) + (santa_map(santa_path) - 1) +  1


print(robo_santa_map('^v'))
print(robo_santa_map('^>v<'))
print(robo_santa_map('^v^v^v^v^v'))
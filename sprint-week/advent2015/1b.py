count = 0
floor = 0

targetString = "(()))"

for char in targetString:
    if floor < 0:
        print(count)
    elif char == "(":
        floor += 1
    elif char == ")":
        floor -= 1
    count += 1

if floor < 0:
    print(count)
else:
    print("never basement")
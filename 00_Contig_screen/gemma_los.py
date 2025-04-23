import argparse

# 参数解析
parser = argparse.ArgumentParser(description="Filter config names based on fixed total length threshold.")
parser.add_argument("input_file", help="Path to the input TSV file")
parser.add_argument("--length", type=int, required=True, help="Threshold length to  aligned printing")

args = parser.parse_args()

input_file = args.input_file
length_threshold = args.length

a = 0
c = 0

with open(input_file) as file_object:
    for line in file_object:
        line = line.strip()
        if not line:
            continue
        line2 = line.split("\t")
        if a > 0:
            if line2[0] == config_name:
                d_value += int(line2[3]) - int(line2[2])
                if c == 0 and d_value > length_threshold:
                    print(config_name)
                    c = 1
            else:
                config_name = line2[0]
                d_value = int(line2[3]) - int(line2[2])
                c = 0
                if d_value > length_threshold:
                    print(config_name)
                    c = 1
        else:
            a += 1
            config_name = line2[0]
            d_value = int(line2[3]) - int(line2[2])
            if d_value > length_threshold:
                print(config_name)
                c = 1

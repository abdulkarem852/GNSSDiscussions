import datetime
from bit_generator import Bitgenerator



rinex_file = "GODS00USA_R_20240830000_01D_GN.rnx"
alm_file = "gpsAlmanac.txt"
time = datetime.datetime(2024, 3, 23, 2, 0, 0, 0)
prn = 5
message = "mohanad is awesome"

bit_generator = Bitgenerator(rinex_file, alm_file, time, prn)

navigation_message = bit_generator.gen_data(message=message)


with open("Navigation Message Bitstream.txt", "w", newline="") as file:
    for frame in range(1, 26):
        for sf in range(1,6):
            file.write(f"F-{frame} SF-{sf} \n")
            for w in range(1, 11):
                ws = (frame - 1) * 1500 + (sf - 1) * 300 + (w-1) * 30
                we = ws + 30
                word = [str(a) for a in navigation_message[ws:we]]
                file.write("".join(word) + '\n')
            file.write("\n")
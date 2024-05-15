import datetime
from bit_generator import Bitgenerator
import argparse

desc = "This program takes a time, a rinex file and a SEM almanac file and uses them to generate a full navigation message bit stream"
parser = argparse.ArgumentParser(description=desc)



parser.add_argument("-r", "--rinex_path", help="Path of the RINEX file", type=str)
parser.add_argument("-a", "--almanac", help="Path of the SEM Almanac file", type=str)
parser.add_argument("-t", "--time", help="YYYY-MM-DDTHH:MM:SS", type=str)
parser.add_argument("-p", "--prn", help="Enter the PRN of the satellite that will transmit the data", type=int)
parser.add_argument("-m", "--message", help="Enter the Message that you'd like to be sent on the Navigation message (22 characters max)(optional)", type=str, default="No Message")
parser.add_argument("-f", "--file_name", help="Name of the file you'd like the navigation message to be stored in (without extension)(default=Navigation Message Bitstream)",type=str , default="Navigation Message Bitstream")

args = parser.parse_args()

time = datetime.datetime.strptime(args.time, "%Y-%m-%dT%H:%M:%S")
rinex_file = args.rinex_path
alm_file = args.almanac
prn = args.prn
message= args.message 
fname = args.file_name + ".txt"

# rinex_file = "GODS00USA_R_20240830000_01D_GN.rnx"
# alm_file = "gpsAlmanac.txt"
# # filename = "brdc0830.24n"
# time = datetime.datetime(2024, 3, 23, 2, 0, 0, 0)
# prn = 5
# # frame = 18
bit_generator = Bitgenerator(rinex_file, alm_file, time, prn)
# message = "mohanad is awesome"

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




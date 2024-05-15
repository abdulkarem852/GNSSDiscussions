import math
import georinex as gr
import datetime
import numpy as np
import pandas as pd
import warnings


GPS_EPOCH = datetime.datetime(1980, 1, 6, 0, 0, 0, 0)
LEAP_SECONDS = 18
SECS_PER_WEEK = (7 * 24 * 3600)
MU = 3.986005E+14
C = 2.99792458E+8
AREF = 26559710
J_ITER = 5
REL_FACTOR = -2 * np.sqrt(MU) / (C*C)
xPI = 3.1415926535898
OmegaDotE = 7.2921151467E-5
OmegaDotREF = -2.6E-9

PREAMBLE = [1,0,0,0,1,0,1,1]
TLM_WORD = PREAMBLE + [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0]

def dec_bin(num, bits, scale_factor=1):
    numx = round(num / scale_factor)
    if numx > 2**bits:
        raise ValueError(f"{numx} is too large to be represented in {bits} bits")
    if numx >= 0:
        binary_conv = bin(numx)[2:].zfill(bits)
    else:
        binary_conv = bin((1 << bits) + numx)[2:]
    
    return [int(n) for n in binary_conv]


def to_ascii(text):
    letters_ascii = [ord(letter) for letter in text]

    letters_ascii_bin = []
    for letter in letters_ascii:
        letters_ascii_bin += dec_bin(letter, bits=8)
    
    if len(letters_ascii_bin) > 176:
        return letters_ascii_bin[:176]
    
    while len(letters_ascii_bin) != 176:
        letters_ascii_bin += dec_bin(32, bits=8)
    
    return letters_ascii_bin



class Bitgenerator:
    def __init__(self, rinex_file,alm_file, time, prn) -> None:
        self.sv = "G%02d" % prn
        self.time = time

        self.rinex_file = rinex_file
        self.eph = self.readRinexFile()
        self.sv_eph = self.eph[self.sv]

        self.alm_file = alm_file
        self.alm = self.readSemAlmanac()
        
        self.how_tow = round(self.sv_eph["TransTime"] / 6)


    def readSemAlmanac(self):
        with open(self.alm_file, "r") as file:
            lines = file.readlines()

        lines = [a.replace("\n", "").strip() for a in lines]

        temp_list = []
        sections = []
        for line in lines:
            if line == "":
                sections.append(temp_list)
                temp_list = []
            else:
                temp_list.append(line)
        sections.append(temp_list)
        almanac = {
            "num_svs"   :   int(sections[0][0][:2]),
            "WNa"       :   int(sections[0][1][:3]),
            "toa"       :   int(sections[0][1][4:]),
        }

        sections = sections[1:]
        for section in sections:
            line1 = [float(x) for x in section[3].split()]
            line2 = [float(x) for x in section[4].split()]
            line3 = [float(x) for x in section[5].split()]
            sv_alm = {
                "id"        :   int(section[0]),
                "SVID"      :   int(section[1]),
                "URA"       :   int(section[2]),
                "e"         :   line1[0],
                "delta_i"   :   line1[1],
                "OmegaDot"  :   line1[2],
                "sqrtA"     :   line2[0],
                "Omega0"    :   line2[1],
                "omega"     :   line2[2],
                "M0"        :   line3[0],
                "Af0"       :   line3[1],
                "Af1"       :   line3[2],
                "health"    :   int(section[6]),
                "config"    :   int(section[7]),   
            } 
            almanac["G%02d" % sv_alm["id"]] = sv_alm
        
        return almanac


    def readRinexFile(self):
        rinex_nav_file = gr.load(self.rinex_file)
        user_tx = pd.Timestamp(datetime.datetime.strftime(self.time, "%Y-%m-%d %H:%M:%S")) 

        eph = {}

        # Load ephemeris for all space vehicles at user time
        for svn in range(1,33):
            sv_prn = "G%02d" % svn
            eph_t = rinex_nav_file.sel(sv=sv_prn, time=user_tx)
            # n=0
            # while math.isnan(eph_t.isel(time=n)["IODE"]):
            #     n+=1
            # eph_t = eph_t.isel(time=n)
            iode = float(eph_t["IODE"]) 
            toc = (user_tx - GPS_EPOCH).total_seconds() % SECS_PER_WEEK + LEAP_SECONDS
            if math.isnan(iode):
                # print("PRN not found: ", sv_prn)
                continue
            else:
                
                keyx = list(eph_t.keys())
                eph_x = {
                    "Af0"           :   float(eph_t["SVclockBias"]),
                    "Af1"           :   float(eph_t["SVclockDrift"]),
                    "Af2"           :   float(eph_t["SVclockDriftRate"]),
                    "toc"           :   toc,
                    "C_rs"          :   float(eph_t[keyx[4]]),
                    "Delta_n"       :   float(eph_t[keyx[5]]),
                    "M0"            :   float(eph_t[keyx[6]]),
                    "C_uc"          :   float(eph_t[keyx[7]]), 
                    "e"             :   float(eph_t[keyx[8]]),
                    "C_us"          :   float(eph_t[keyx[9]]),
                    "sqrtA"         :   float(eph_t[keyx[10]]),
                    "toe"           :   float(eph_t["Toe"]),
                    "C_ic"          :   float(eph_t[keyx[12]]),
                    "Omega"         :   float(eph_t[keyx[13]]),
                    "C_is"          :   float(eph_t[keyx[14]]),
                    "i0"            :   float(eph_t[keyx[15]]),
                    "C_rc"          :   float(eph_t[keyx[16]]),
                    "omega"         :   float(eph_t[keyx[17]]),
                    "OmegaDot"      :   float(eph_t[keyx[18]]),
                    "IDOT"          :   float(eph_t[keyx[19]]),
                    "GPSWeek"       :   float(eph_t[keyx[21]]),
                    "IODE"          :   iode,
                    "IODC"          :   float(eph_t["IODC"]),
                    "L2PDataFlag"   :   float(eph_t["L2Pflag"]),
                    "L2ChannelCode" :   float(eph_t["CodesL2"]),
                    "SVAcc"         :   float(eph_t["SVacc"]),
                    "health"        :   float(eph_t["health"]),
                    "TGD"           :   float(eph_t["TGD"]),
                    "GPSWeek"       :   float(eph_t["GPSWeek"]),
                    "FitInterval"   :   float(eph_t["FitIntvl"]),
                    "TransTime"     :   float(eph_t["TransTime"]),
                }
                eph[sv_prn] = eph_x
        
        return eph
    

    def gen_word(self, bits, D30_star, D29_star):
        parity_bits = self.gen_p(bits, D29_star=D29_star, D30_star=D30_star)

        if D30_star == 1:
            return (np.bitwise_not(bits) + 2).tolist() + parity_bits

        return bits + parity_bits


    def gen_t(self, bits1_22, D29_star, D30_star):
        # if D30_star == 1:
        #     bits1_22 = np.bitwise_not(bits1_22) + 2
        tmp_D29 = D30_star ^ bits1_22[0] ^ bits1_22[2] ^ bits1_22[4] ^ bits1_22[5] ^ bits1_22[6] ^ bits1_22[8] ^ bits1_22[9] ^ bits1_22[13] ^ bits1_22[14] ^ bits1_22[15] ^ bits1_22[16] ^ bits1_22[17] ^ bits1_22[20] ^ bits1_22[21]  
        d24 = tmp_D29
        tmp_D30 = D29_star ^ bits1_22[2] ^ bits1_22[4] ^ bits1_22[5] ^ bits1_22[7] ^ bits1_22[8] ^ bits1_22[9] ^ bits1_22[10] ^ bits1_22[12] ^ bits1_22[14] ^ bits1_22[18] ^ bits1_22[21] ^ d24
        d23 = tmp_D30
        return [d23, d24] 


    def gen_p(self, bits_1_24, D29_star, D30_star):
        # if D30_star == 1:
        #     bits_1_24 = np.bitwise_not(bits_1_24)
        #     bits_1_24 += 2
        D25 = D29_star ^ bits_1_24[0] ^ bits_1_24[1] ^ bits_1_24[2] ^ bits_1_24[4] ^ bits_1_24[5] ^ bits_1_24[9] ^ bits_1_24[10] ^ bits_1_24[11] ^ bits_1_24[12] ^ bits_1_24[13] ^ bits_1_24[16] ^ bits_1_24[17] ^ bits_1_24[19] ^ bits_1_24[22]
        D26 = D30_star ^ bits_1_24[1] ^ bits_1_24[2] ^ bits_1_24[3] ^ bits_1_24[5] ^ bits_1_24[6] ^ bits_1_24[10] ^ bits_1_24[11] ^ bits_1_24[12] ^ bits_1_24[13] ^ bits_1_24[14] ^ bits_1_24[17] ^ bits_1_24[18] ^ bits_1_24[20] ^ bits_1_24[23]
        D27 = D29_star ^ bits_1_24[0] ^ bits_1_24[2] ^ bits_1_24[3] ^ bits_1_24[4] ^ bits_1_24[6] ^ bits_1_24[7] ^ bits_1_24[11] ^ bits_1_24[12] ^ bits_1_24[13] ^ bits_1_24[14] ^ bits_1_24[15] ^ bits_1_24[18] ^ bits_1_24[19] ^ bits_1_24[21]
        D28 = D30_star ^ bits_1_24[1] ^ bits_1_24[3] ^ bits_1_24[4] ^ bits_1_24[5] ^ bits_1_24[7] ^ bits_1_24[8] ^ bits_1_24[12] ^ bits_1_24[13] ^ bits_1_24[14] ^ bits_1_24[15] ^ bits_1_24[16] ^ bits_1_24[19] ^ bits_1_24[20] ^ bits_1_24[22]
        D29 = D30_star ^ bits_1_24[0] ^ bits_1_24[2] ^ bits_1_24[4] ^ bits_1_24[5] ^ bits_1_24[6] ^ bits_1_24[8] ^ bits_1_24[9] ^ bits_1_24[13] ^ bits_1_24[14] ^ bits_1_24[15] ^ bits_1_24[16] ^ bits_1_24[17] ^ bits_1_24[20] ^ bits_1_24[21] ^ bits_1_24[23]
        D30 = D29_star ^ bits_1_24[2] ^ bits_1_24[4] ^ bits_1_24[5] ^ bits_1_24[7] ^ bits_1_24[8] ^ bits_1_24[9] ^ bits_1_24[10] ^ bits_1_24[12] ^ bits_1_24[14] ^ bits_1_24[18] ^ bits_1_24[21] ^ bits_1_24[22] ^ bits_1_24[23]

        return [D25, D26, D27, D28, D29, D30]


    def gen_how(self, subframe, frame):
        how_tow = self.how_tow + (frame - 1) * 5 + subframe 
        how_tow_bin = dec_bin(how_tow, bits=17)
        alert_flag = [0]
        as_flag = [0]
        sf_id = dec_bin(subframe, bits=3)
        bits_1_22 = how_tow_bin + alert_flag + as_flag + sf_id
        bits_23_24 = self.gen_t(bits1_22=bits_1_22, D29_star=TLM_WORD[-2], D30_star=TLM_WORD[-1])
        bits_1_24 = bits_1_22 + bits_23_24
        how = self.gen_word(bits_1_24, D29_star=TLM_WORD[-2], D30_star=TLM_WORD[-1])
        return how
    

    def gen_sv_health(self, prn):

        sv = "G%02d" % prn
        sv_health = self.alm[sv]["health"] if sv in self.alm else 63
        return dec_bin(sv_health, bits=6)


    def gen_sv_config(self, prn):
        sv = "G%02d" % prn
        sv_config = self.alm[sv]["config"] if sv in self.alm else 15
        return dec_bin(sv_config, bits=4)


    def gen_alm_sv(self, prn, how):
        data_id = [0, 1]
        sv = "G%02d" % prn

        if sv not in self.alm:
            sv_id = [0] * 6
            dummy_data_w3 = self.gen_word(data_id + sv_id +[1,0] * 8, D29_star=how[-2], D30_star=how[-1])
            dummy_data_w4_9 = self.gen_word([1,0]*12, D29_star=dummy_data_w3[-2], D30_star=dummy_data_w3[-1]) * 6
            dummy_data_w10 = self.gen_word([1,0]*11 + self.gen_t([1,0]*11, D29_star=dummy_data_w4_9[-2], D30_star=dummy_data_w4_9[-1]), D29_star=dummy_data_w4_9[-2], D30_star=dummy_data_w4_9[-1])
            dummy_data = dummy_data_w3 + dummy_data_w4_9 + dummy_data_w10
            warnings.warn(f"{sv} isn't available in the almanac. Dummy data was filled in its place") 
            return TLM_WORD + how + dummy_data 
        
        sv_id = dec_bin(prn, bits=6)
        # sv_id = dec_bin(self.alm[sv]["SVID"], bits=6)

        e = dec_bin(self.alm[sv]["e"], bits=16, scale_factor=2**-21)
        toa = dec_bin(self.alm["toa"], bits=8, scale_factor=2**12)
        di = dec_bin(self.alm[sv]["delta_i"] , bits=16, scale_factor=2**-19)
        OmegaDot = dec_bin(self.alm[sv]["OmegaDot"] , bits=16, scale_factor=2**-38)
        sv_health = dec_bin(self.alm[sv]["health"], bits=8)
        sv_health = [0] * 3 + sv_health[3:]
        sqrtA = dec_bin(self.alm[sv]["sqrtA"], bits=24, scale_factor=2**-11)
        Omega0 = dec_bin(self.alm[sv]["Omega0"] , bits=24, scale_factor=2**-23)
        omega= dec_bin(self.alm[sv]["omega"] , bits=24, scale_factor=2**-23)
        m0 = dec_bin(self.alm[sv]["M0"] , bits=24, scale_factor=2**-23)
        af0 = dec_bin(self.alm[sv]["Af0"], bits=11, scale_factor=2**-20)
        af1 = dec_bin(self.alm[sv]["Af1"], bits=11, scale_factor=2**-38)

        w3 = self.gen_word(data_id + sv_id + e, D29_star=how[-2], D30_star=how[-1])
        w4 = self.gen_word(toa + di, D29_star=w3[-2], D30_star=w3[-1])
        w5 = self.gen_word(OmegaDot + sv_health, D29_star=w4[-2], D30_star=w4[-1])
        w6 = self.gen_word(sqrtA, D29_star=w5[-2], D30_star=w5[-1])
        w7 = self.gen_word(Omega0, D29_star=w6[-2], D30_star=w6[-1])
        w8 = self.gen_word(omega, D29_star=w7[-2], D30_star=w7[-1])
        w9 = self.gen_word(m0, D29_star=w8[-2], D30_star=w8[-1])
        
        w10 = self.gen_word(af0[:8] + af1 + af0[8:] + self.gen_t(af0[:8] + af1 + af0[8:], D29_star=w9[-2], D30_star=w9[-1]), D29_star=w9[-2], D30_star=w9[-1])

        return TLM_WORD + how + w3 + w4 + w5 + w6 + w7 + w8 + w9 + w10
    

    def gen_sf1(self, frame):                                   
        how = self.gen_how(subframe=1, frame=frame)

        gps_week_bin = dec_bin(int(self.sv_eph["GPSWeek"]) % 1024, bits=10)
        w3_bits_11_22 = dec_bin(self.sv_eph["L2ChannelCode"], bits=2) + [0] * 10
        IODC = dec_bin(self.sv_eph["IODC"], bits=10)
        w3_bits_1_24 = gps_week_bin + w3_bits_11_22 + IODC[:2]
        w3 = self.gen_word(w3_bits_1_24, D29_star=how[-2], D30_star=how[-1])

        w4_bits_1_24 = [int(self.sv_eph["L2PDataFlag"])] + [0] * 23
        w4 = self.gen_word(bits=w4_bits_1_24, D29_star=w3[-2], D30_star= w3[-1])

        w5 = self.gen_word([0] * 24, D29_star=w4[-2], D30_star=w4[-1])

        w6 = self.gen_word([0] * 24, D29_star=w5[-2], D30_star=w5[-1])

        w7_bits_1_24 = [0] * 16 + dec_bin(self.sv_eph["TGD"], bits=8, scale_factor=2**-31)
        w7 = self.gen_word(bits=w7_bits_1_24, D29_star=w6[-2], D30_star= w6[-1])

        w8_bits_1_24 = IODC[2:] + dec_bin(self.sv_eph["toc"], bits=16, scale_factor= 2**4)
        w8 = self.gen_word(bits=w8_bits_1_24, D29_star=w7[-2], D30_star= w7[-1])

        af2_bin = dec_bin(self.sv_eph["Af2"], bits=8,scale_factor=2**-55)
        af1_bin = dec_bin(self.sv_eph["Af1"], bits=16, scale_factor=2**-43)
        af0_bin = dec_bin(self.sv_eph["Af0"], bits=22, scale_factor=2**-31)
        w9 = self.gen_word(af2_bin + af1_bin, D29_star=w8[-2], D30_star=w8[-1])

        w10_bits_1_24 = af0_bin + self.gen_t(af0_bin, D29_star=w9[-2], D30_star=w9[-1])
        w10 = self.gen_word(bits=w10_bits_1_24, D29_star=w9[-2], D30_star= w9[-1])

        return TLM_WORD + how + w3 + w4 + w5 + w6 + w7 + w8 + w9 + w10    


    def gen_sf2(self, frame):
        how = self.gen_how(subframe=2, frame=frame)

        w3_bits_1_24 = dec_bin(self.sv_eph["IODE"], bits=8) + dec_bin(self.sv_eph["C_rs"], bits=16, scale_factor=2**-5)
        w3 = self.gen_word(bits= w3_bits_1_24, D29_star=how[-2], D30_star=how[-1])

        m0 = dec_bin(self.sv_eph["M0"] / xPI, bits=32, scale_factor=2**-31)
        w4_bits_1_24 = dec_bin(self.sv_eph["Delta_n"] / xPI, bits= 16, scale_factor=2**-43) + m0[:8]
        w4 = self.gen_word(w4_bits_1_24, D29_star=w3[-2], D30_star=w3[-1])

        w5_bits_1_24 = m0[8:]
        w5 = self.gen_word(w5_bits_1_24, D29_star=w4[-2], D30_star=w4[-1])

        e = dec_bin(self.sv_eph["e"], bits=32, scale_factor=2**-33)
        w6_bits_1_24 = dec_bin(self.sv_eph["C_uc"], bits=16, scale_factor=2**-29) + e[:8]
        w6 = self.gen_word(w6_bits_1_24, D29_star= w5[-2], D30_star=w5[-1])

        w7_bits_1_24 = e[8:]
        w7 = self.gen_word(w7_bits_1_24, D29_star=w6[-2], D30_star=w6[-1])

        sqrtA = dec_bin(self.sv_eph["sqrtA"], bits=32, scale_factor=2**-19)
        w8_bits_1_24 = dec_bin(self.sv_eph["C_us"], bits=16, scale_factor=2**-29) + sqrtA[:8]
        w8 = self.gen_word(w8_bits_1_24, D29_star=w7[-2], D30_star=w7[-1])

        w9_bits_1_24 = sqrtA[8:]
        w9 = self.gen_word(w9_bits_1_24, D29_star=w8[-2], D30_star=w8[-1])

        w10_bits_1_16 = dec_bin(self.sv_eph["toe"], bits=16, scale_factor=2**4)
        fit_flag = [1] if self.sv_eph["FitInterval"] > 4 else [0]
        AODO = 0 #TODO: Look for what it is upposed to be, 0 is just a place holder.
        AODO_bin = dec_bin(AODO, bits=5, scale_factor=900)
        w10_bits_1_22 = w10_bits_1_16 + fit_flag + AODO_bin
        w10_bits_1_24 = w10_bits_1_22 + self.gen_t(w10_bits_1_22, D29_star=w9[-2], D30_star=w9[-1])
        w10 = self.gen_word(w10_bits_1_24, D29_star=w9[-2], D30_star=w9[-1])

        return TLM_WORD + how + w3 + w4 + w5 + w6 + w7 + w8 + w9 + w10   
    

    def gen_sf3(self, frame):
        how = self.gen_how(subframe=3, frame=frame)
        
        Omega0 = dec_bin(self.sv_eph["Omega"] / xPI, bits=32, scale_factor=2**-31)
        w3_bits_1_24 = dec_bin(self.sv_eph["C_ic"], bits=16, scale_factor=2**-29) + Omega0[:8]
        w3 = self.gen_word(w3_bits_1_24, D29_star=how[-2], D30_star=how[-1])

        w4_bits_1_24 = Omega0[8:]
        w4 = self.gen_word(w4_bits_1_24, D29_star=w3[-2], D30_star=w3[-1])

        i0 = dec_bin(self.sv_eph["i0"] / xPI, bits=32, scale_factor=2**-31)
        w5_bits_1_24 = dec_bin(self.sv_eph["C_is"], bits=16, scale_factor=2**-29) + i0[:8]
        w5 = self.gen_word(w5_bits_1_24, D29_star=w4[-2], D30_star=w4[-1])

        w6_bits_1_24 = i0[8:]
        w6 = self.gen_word(w6_bits_1_24, D29_star=w5[-2], D30_star=w5[-1])

        omega = dec_bin(self.sv_eph["omega"] / xPI, bits=32, scale_factor=2**-31)
        w7_bits_1_24 = dec_bin(self.sv_eph["C_rc"], bits=16, scale_factor=2**-5) + omega[:8]
        w7 = self.gen_word(w7_bits_1_24, D29_star=w6[-2], D30_star=w6[-1])

        w8_bits_1_24 = omega[8:]
        w8 = self.gen_word(w8_bits_1_24, D29_star=w7[-2], D30_star=w7[-1])

        w9_bits_1_24 = dec_bin(self.sv_eph["OmegaDot"] / xPI, bits=24, scale_factor=2**-43)
        w9 = self.gen_word(w9_bits_1_24, D29_star=w8[-2], D30_star=w8[-1])

        w10_bits_1_22 = dec_bin(self.sv_eph["IODE"], bits=8) + dec_bin(self.sv_eph["IDOT"] / xPI, bits=14, scale_factor=2**-43) 
        w10_bits_1_24 = w10_bits_1_22 + self.gen_t(w10_bits_1_22, D29_star=w9[-2], D30_star=w9[-1])
        w10 = self.gen_word(w10_bits_1_24, D29_star=w9[-2], D30_star=w9[-1])

        return TLM_WORD + how + w3 + w4 + w5 + w6 + w7 + w8 + w9 + w10   


    def gen_sf4(self, frame, message="No message sent"):
        data_id = [0,1]                             ### Used to indicate data structure, which is LNAV here. (per my understanding of pg.113 IS-GPS-200 August 2022)
        reserved_pages = {                          ### Refer to pg. 113 IS-GPS-200,, August 2022
            1   :   57,
            6   :   57,
            11  :   57,
            16  :   57,
            21  :   57,
            14  :   53,
            15  :   54,
            12  :   62,
            19  :   58,
            20  :   59,
            22  :   60,
            23  :   61,
            24  :   62,
        }

        how = self.gen_how(subframe=4, frame=frame)

        if frame in reserved_pages:
            sv_id = dec_bin(reserved_pages[frame], bits=6)
            w3 = self.gen_word( data_id + sv_id + [0] * 16, D29_star=how[-2], D30_star=how[-1])

            w4_9 = self.gen_word([0] * 24, D29_star=w3[-2], D30_star=w3[-1])
            for i in range(5):
                w4_9 += self.gen_word([0] * 24, D29_star=w4_9[-2], D30_star=w4_9[-1])
            
            w10_bits_1_24 =[0] * 22 + self.gen_t([0] * 22, D29_star=w4_9[-2], D30_star=w4_9[-1])
            w10 = self.gen_word(w10_bits_1_24, D29_star=w4_9[-2], D30_star=w4_9[-1])

            return TLM_WORD + how + w3 + w4_9 + w10
        
        elif frame == 13:                ### This frame contains the NMCT(Navigation Message Correction Table)
            sv_id = dec_bin(52, bits=6)
            av_indicator = [0, 0]        ### These 2 bits are the availability indicators for the NMCT. Since I couldn't find it in a Rinex file, I will set it to all zeros. See pg. 123 IS-GPS-200, August 2022

            w3 = self.gen_word( data_id + sv_id + av_indicator + [0] * 14, D29_star=how[-2], D30_star=how[-1])

            w4_9 = self.gen_word([0] * 24, D29_star=w3[-2], D30_star=w3[-1])
            for i in range(5):
                w4_9 += self.gen_word([0] * 24, D29_star=w4_9[-2], D30_star=w4_9[-1])
            
            w10_bits_1_24 =[0] * 22 + self.gen_t([0] * 22, D29_star=w4_9[-2], D30_star=w4_9[-1])
            w10 = self.gen_word(w10_bits_1_24, D29_star=w4_9[-2], D30_star=w4_9[-1])

            return TLM_WORD + how + w3 + w4_9 + w10
        
        elif frame == 17:           ### This frame contains a special messagede coded in ASCII
            sv_id = dec_bin(55, bits = 6)

            message = to_ascii(message)

            w3_9 = self.gen_word(data_id + sv_id + message[:16], D29_star=how[-2], D30_star=how[-1])

            for i in range(6):
                w3_9 += self.gen_word(message[(24 * i) + 16 : (24 * (i + 1) + 16)], D29_star=w3_9[-2], D30_star=w3_9[-1])

            w10_bits_1_22 = message[160:] + [0] * 6
            w10_bits_1_24 = w10_bits_1_22 + self.gen_t(w10_bits_1_22, D29_star=w3_9[-2], D30_star=w3_9[-1])
            w10 = self.gen_word(w10_bits_1_24, D29_star=w3_9[-2], D30_star=w3_9[-1])

            return TLM_WORD + how + w3_9 + w10
        
        elif frame == 18:           ### This frame contains Iono parameters, they're not available in every rinex navigation file.
            ############################ Description of these parameters is given in pg.125-128 IS-GPS-200, August 2022 #####################################
            alpha0 = alpha1 = alpha2 = alpha3 = [0] * 8     ### Not provided in current RINEX file we are using, 0 used
            beta0 = beta1 = beta2 = beta3 = [0] * 8         ### Not provided in current RINEX file we are using, 0 used
            A0 = 0                                          ### Not provided in current RINEX file we are using, 0 used
            A1 = 0                                          ### Not provided in current RINEX file we are using, 0 used
            Dt_LS = 18                                      ### Not provided in current RINEX file we are using, 0 used
            tot = 0                                         ### Not provided in current RINEX file we are using, 0 used
            WNt = 2149 % 256              ### Not sure about this
            WN_LSF = WNt % 256                              ### Not sure about this
            DN = 1                                          ### Not sure about this
            Dt_LSF = 18                                     ### Not sure about this
            ############################# Word Generation Algorithm ###########################################################################################
            sv_id = dec_bin(56, bits=6)
            w3 = self.gen_word( data_id + sv_id + alpha0 + alpha1, D29_star=how[-2], D30_star=how[-1])

            w4 = self.gen_word(alpha2 + alpha3 + beta0, D29_star=w3[-2], D30_star=w3[-1])

            w5 = self.gen_word(beta1 + beta2 + beta3, D29_star=w4[-2], D30_star=w4[-1])

            w6 = self.gen_word(dec_bin(A1, bits=24, scale_factor=2**-50), D29_star=w5[-2], D30_star=w5[-1])

            A0_bin = dec_bin(A0, bits=32, scale_factor= 2**-30)
            w7 = self.gen_word(A0_bin[:24], D29_star=w6[-2], D30_star=w6[-1])

            w8 = self.gen_word(A0_bin[24:] + dec_bin(tot, bits=8, scale_factor=2**12) + dec_bin(WNt, bits = 8), D29_star=w7[-2], D30_star=w7[-1])

            w9 = self.gen_word(dec_bin(Dt_LS, bits=8) + dec_bin(WN_LSF,bits=8) + dec_bin(DN, bits=8), D29_star=w8[-2], D30_star=w8[-1])

            w10_bits_1_22 = dec_bin(Dt_LSF, bits=8) + [0] * 14
            w10 = self.gen_word(w10_bits_1_22 + self.gen_t(w10_bits_1_22, D29_star=w9[-2], D30_star=w9[-1]), D29_star=w9[-2], D30_star=w9[-1])

            return TLM_WORD + how + w3 + w4 + w5 + w6 + w7 + w8 + w9 + w10
        
        elif frame == 25:
            sv_id = dec_bin(63, bits=6)
            w3_bits_1_24 = data_id + sv_id
            for i in range(1, 5):
                as_config = self.gen_sv_config(i)
                w3_bits_1_24 += as_config
            
            w3_8 = self.gen_word(w3_bits_1_24, D29_star=how[-2], D30_star=how[-1])

            for i in range(5):
                bits_1_24 = []
                for n in range(6):
                    prn = 5 + n + i * 6 
                    if prn <= 32:
                        as_config = self.gen_sv_config(prn)
                    else:
                        sys_bits = [0] * 2
                        bits_1_24 += sys_bits + self.gen_sv_health(25)
                        break
                    bits_1_24 += as_config
                w3_8 += self.gen_word(bits_1_24, D29_star=w3_8[-2], D30_star=w3_8[-1])
            
            sv_healths = []
            for i in range(7):
                sv_healths += self.gen_sv_health(26 + i)
            
            w9 = self.gen_word(sv_healths[:24], D29_star=w3_8[-2], D30_star=w3_8[-1])

            w10_bits_1_22 = sv_healths[24:] + [0] * 4

            w10_bits_1_24 = w10_bits_1_22 + self.gen_t(w10_bits_1_22, D29_star=w9[-2], D30_star=w9[-1])
            w10 = self.gen_word(w10_bits_1_24, D29_star=w9[-2], D30_star=w9[-1])

            return TLM_WORD + how + w3_8 + w9 + w10
        else:
            prn = frame + 23 if frame < 6 else frame + 22

            return self.gen_alm_sv(prn, how)

    
    def gen_sf5(self, frame):
        how = self.gen_how(subframe=5, frame=frame)

        data_id = [0,1]
        sv_id = dec_bin(51, bits=6)         ### Used only for Frame 25

        if frame < 25:
            return self.gen_alm_sv(frame, how)
        else:
            w3 = self.gen_word(data_id + sv_id + dec_bin(self.alm["toa"], bits=8, scale_factor=2**12) + dec_bin(self.alm["WNa"], bits=8), D29_star=how[-2], D30_star=how[-1])

            sv_healths = []
            for i in range(1, 25):
                sv_healths += self.gen_sv_health(i)
            
            w3_9 = w3
            for i in range(6):
                w3_9 += self.gen_word(sv_healths[i * 24:(i + 1) * 24], D29_star=w3_9[-2], D30_star=w3_9[-1])
            
            w10 = self.gen_word([0] * 22 + self.gen_t([0] * 22, D29_star=w3_9[-2], D30_star=w3_9[-1]), D29_star=w3_9[-2], D30_star=w3_9[-1])
            
            return TLM_WORD + how + w3_9 + w10
            

    def gen_frame(self, frame, message="No message sent"):
        return self.gen_sf1(frame) + self.gen_sf2(frame) + self.gen_sf3(frame) + self.gen_sf4(frame, message) + self.gen_sf5(frame) 
    

    def gen_data(self, message="No message sent"):
        data = []
        for i in range(1, 26):
            data += self.gen_frame(i, message=message)
        return data

    



            



        



        
            

        


            




            


        

            




            
    












                                  
    




        




with open("sample_input_zz_l04mu.txt", "w", encoding="utf-8") as file:
    for i in range(0,390):
        if i<10: file.write("/cefs/higgs/liugeliang/CEPC/202501/Production/4fermions/zz_l04mu/Combined/rec_zz_l04mu_0000"+str(i)+".root\n")
        else: file.write("/cefs/higgs/liugeliang/CEPC/202501/Production/4fermions/zz_l04mu/Combined/rec_zz_l04mu_000"+str(i)+".root\n")

with open("sample_input_zz_l0mumu.txt", "w", encoding="utf-8") as file:
    for i in range(0,390):
        if i<10: file.write("/cefs/higgs/liugeliang/CEPC/202501/Production/4fermions/zz_l0mumu/Combined/rec_zz_l0mumu_0000"+str(i)+".root\n")
        else: file.write("/cefs/higgs/liugeliang/CEPC/202501/Production/4fermions/zz_l0mumu/Combined/rec_zz_l0mumu_000"+str(i)+".root\n")

with open("sample_input_zzorww_l0mumu.txt", "w", encoding="utf-8") as file:
    for i in range(0,390):
        if i<10: file.write("/cefs/higgs/liugeliang/CEPC/202501/Production/4fermions/zzorww_l0mumu/Combined/rec_zzorww_l0mumu_0000"+str(i)+".root\n")
        else: file.write("/cefs/higgs/liugeliang/CEPC/202501/Production/4fermions/zzorww_l0mumu/Combined/rec_zzorww_l0mumu_000"+str(i)+".root\n")


with open("sample_input_zz_sl0mu_down.txt", "w", encoding="utf-8") as file:
    for i in range(0,390):
        if i<10: file.write("/cefs/higgs/liugeliang/CEPC/202501/Production/4fermions/zz_sl0mu_down/Combined/rec_zz_sl0mu_down_0000"+str(i)+".root\n")
        else: file.write("/cefs/higgs/liugeliang/CEPC/202501/Production/4fermions/zz_sl0mu_down/Combined/rec_zz_sl0mu_down_000"+str(i)+".root\n")

with open("sample_input_zz_sl0mu_up.txt", "w", encoding="utf-8") as file:
    for i in range(0,390):
        if i<10: file.write("/cefs/higgs/liugeliang/CEPC/202501/Production/4fermions/zz_sl0mu_up/Combined/rec_zz_sl0mu_up_0000"+str(i)+".root\n")
        else: file.write("/cefs/higgs/liugeliang/CEPC/202501/Production/4fermions/zz_sl0mu_up/Combined/rec_zz_sl0mu_up_000"+str(i)+".root\n")

with open("sample_input_ll.txt", "w", encoding="utf-8") as file:
    for i in range(0,974):
        if i<10: file.write("/cefs/higgs/zhangkl/Production/2410/E240_2f_ll/Combined/rec_E240_2f_ll_0000"+str(i)+".root\n")
        elif i<100: file.write("/cefs/higgs/zhangkl/Production/2410/E240_2f_ll/Combined/rec_E240_2f_ll_000"+str(i)+".root\n")
        else: file.write("/cefs/higgs/zhangkl/Production/2410/E240_2f_ll/Combined/rec_E240_2f_ll_00"+str(i)+".root\n")
with open("sample_input_sznu_l0mumu.txt", "w", encoding="utf-8") as file:
    for i in range(0,99):
        if i<10: file.write("/cefs/higgs/liugeliang/CEPC/202501/Production/4fermions/sznu_l0mumu/Combined/rec_sznu_l0mumu_0000"+str(i)+".root\n")
        else: file.write("/cefs/higgs/liugeliang/CEPC/202501/Production/4fermions/sznu_l0mumu/Combined/rec_sznu_l0mumu_000"+str(i)+".root\n")
with open("sample_input_sig.txt", "w", encoding="utf-8") as file:
    for i in range(0,200):
        if i<10: file.write("/cefs/higgs/zhangkl/Production/2410/E240_mmHinclusive/Combined/rec_E240_mmHinclusive_0000"+str(i)+".root\n")
        elif i<100: file.write("/cefs/higgs/zhangkl/Production/2410/E240_mmHinclusive/Combined/rec_E240_mmHinclusive_000"+str(i)+".root\n")
        else: file.write("/cefs/higgs/zhangkl/Production/2410/E240_mmHinclusive/Combined/rec_E240_mmHinclusive_00"+str(i)+".root\n")

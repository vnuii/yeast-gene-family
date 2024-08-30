# -*- coding: utf-8 -*-
import re
import sys
import os
import subprocess

# change your path to cafe5
cafe5_path = "/PATH/TO/CAFE5"

# os.environ['OPENBLAS_NUM_THREADS'] = '1'
o_filename = []
p_filename = []
treename = []
p_final_input={}

def cafe5_auto(o,p,t,output_file_path,sp,hold_tree_file):
	#make tree
	get_tree_file(t,sp,hold_tree_file)

	#get file path
	o_filename = sorted(os.listdir(o))
	p_filename = sorted(os.listdir(p))
	treename = sorted(os.listdir(t))
	
	#final_output = f"{workdir}/final.o"
	#if not os.path.exists(final_output):
	#    os.mkdir(final_output)
	#final_output_noz = f"{workdir}/final_noz.o"
	#if not os.path.exists(final_output_noz):
	#    os.mkdir(final_output_noz)
	new_final_output = f"{workdir}/new_final.o"
	if not os.path.exists(new_final_output):
		os.mkdir(new_final_output)
	forAL_path = f"{workdir}/forAL_path"
	if not os.path.exists(forAL_path):
		os.mkdir(forAL_path)

	if len(p_filename)!= len(treename):
		print("Error: number of file is not same!")
	else:
		for file1, file2, file3 in zip(p_filename, treename,o_filename):
			pog_input=f"{p}/{file1}"
			t_input = f"{t}/{file2}"
			f2= file2.split('.')
			err_output_path = errormodel(pog_input,t_input,f2)
			k_path = k_parament(pog_input,t_input,f2)
			k_output_path = k_path[0]
			k_tmpdir_path = k_path[1]
			og_input = f"{o}/{file3}"
			# 指定传入R里的参数
			kout_dir_path = k_output_path
			OGcounts_filtered_path = k_tmpdir_path # 删过基因家族的表
			#print(OGcounts_filtered_path)
			tree_path = t_input
			erroModel_path = err_output_path
			forAL_path_1 = f"{forAL_path}/{f2[0]}" # 输出文件位置
			OGcounts_path = og_input
			#final_outPath = f"{final_output}/{f2[0]}" # 输出文件位置
			command = cafe5_path
			parament_file_path = pog_input
			#final_outPath_noz = f"{final_output_noz}/{f2[0]}"
			new_final_outpath = f"{new_final_output}/{f2[0]}"
			subprocess.call(['Rscript', 'automate_CAFE_out.r',
				 kout_dir_path, OGcounts_filtered_path, tree_path, erroModel_path,
				 forAL_path_1, OGcounts_path, new_final_outpath, command, parament_file_path])
	
#run errormodel
def errormodel(pog_input,t_input,f2):
	err_pog_input = pog_input
	err_del_og = {}
	str = "Done!"
	for i in range(1,21):
		err_result= []
		err_command = f"{cafe5_path} -i {err_pog_input} -t {t_input} -c 80 -p -e -o {err_model}/{f2[0]}/ >> {log}/{f2[0]}_errmodel.log 2>&1"
		if not os.path.exists(f"{log}/{f2[0]}_errmodel.log"):
			os.system(err_command)
		with open(f"{log}/{f2[0]}_errmodel.log", "r") as errmodel_log, open(f"{err_pog_input}", "r") as p_input:
			success = 1
			for line in errmodel_log.readlines():
				l1 = line.split(':')
				og = re.search('OG[0-9]+', l1[0])
				if og != None:
					tog = og.group()
					err_del_og[tog] = tog
				if str in line:
					print('cafe error model done!')
					errmodel_log.close()
					p_input.close()
					success = 0
			if success == 0:
				break
			else:
				os.system(err_command)
			for line in p_input.readlines():
				if line.startswith('Desc'):
					line = line.split('\t')
					err_result.append(line)
					continue
				line = line.split('\t')
				OG = re.search('OG[0-9]+', line[0]).group()
				if OG != None:
					if OG not in err_del_og:
						err_result.append(line)
		err_pog_input = f"{tmp}/{f2[0]}.err.tmp"
		with open(err_pog_input, "w") as tmpfile:
			for item in err_result:
				tmpfile.write('\t'.join(item))
	err_path = f"{err_model}/{f2[0]}"
	return(err_path)

def k_parament(pog_input,t_input,f2):
	#runing k parament
	k_pog_input = pog_input
	k_del_og = {}
	k_output = f"{k_o}/{f2[0]}"
	if not os.path.exists(k_output):
		os.mkdir(k_output)
	k_tmp = f"{k_t}/{f2[0]}"
	if not os.path.exists(k_tmp):
		os.mkdir(k_tmp)
	str = "Done!"
	for i in range(2,11):
		for t in range(1,21):
			k_result = []
			k_command = f"{cafe5_path} -i {k_pog_input} -t {t_input} -c 80 -p -e{err_model}/{f2[0]}/Base_error_model.txt -k {i} -o {k_output}/{i} >> {k_p_log}/{f2[0]}.k{i}.log 2>&1"
			if not os.path.exists(f"{k_p_log}/{f2[0]}.k{i}.log"):
				os.system(k_command)
			with open(f"{k_p_log}/{f2[0]}.k{i}.log", "r") as k_log, open(f"{k_pog_input}", "r") as p_input:
				success = 1
				for line in k_log.readlines():
					l1 = line.split(':')
					og = re.search('OG[0-9]+', l1[0])
					if og != None:
						tog = og.group()
						k_del_og[tog] = tog
					if str in line:
						print('cafe para-k done!')
						k_log.close()
						p_input.close()
						success = 0
				if success == 0:
					break
				else:
					os.system(k_command)
				for line in p_input.readlines():
					if line.startswith('Desc'):
						line = line.split('\t')
						k_result.append(line)
						continue
					line = line.split('\t')
					OG = re.search('OG[0-9]+', line[0]).group()
					if OG != None:
						if OG not in k_del_og:
							k_result.append(line)
			k_pog_input = f"{k_tmp}/k{i}.tmp"
			with open(f"{k_tmp}/k{i}.tmp", "w") as tmpfile:
				for item in k_result:
					tmpfile.write('\t'.join(item))
	k_path = f"{k_output}"
	k_tmp_path = f"{k_tmp}"
	#print('k_tmp_path:',k_tmp_path)
	return k_path,k_tmp_path

def get_tree_file(t,sp,htf):
	sp_name = []
	sp_name = os.listdir(sp)
    for f in sp_name :
		new_file = f.split('.')
		command = f"gotree prune -i {htf} -o {t}/{new_file[0]}.nh -r `cat {sp}/{f}`"
		os.system(command)
		print("success to make a tree")



if __name__ == '__main__':
	if len(sys.argv) != 7:
		print('usage: python cafe_atuo_final.py origin_og_table_dir parament_og_table_dir node_tree_file_dir output_file_path sp_name hold_tree_file')
		sys.exit(1)

	o = sys.argv[1]
	p = sys.argv[2]
	t = sys.argv[3]
	output_file_path = sys.argv[4]
	sp_name = sys.argv[5]
	hold_tree_file = sys.argv[6]
	
	od = os.path.abspath(o)
	if not os.path.exists(od):
		print('input o_file_path not found.')
		sys.exit(1)
	pd = os.path.abspath(p)
	if not os.path.exists(pd):
		print('input p_file_path not found.')
		sys.exit(1)
	td = os.path.abspath(t)
	if not os.path.exists(td):
		os.mkdir(td)
		print('input tree_file_path not found.')
		#sys.exit(1)
	if not os.path.exists(sp_name):
		print('input sp_name not found.')
		sys.exit(1)
	if not os.path.isfile(hold_tree_file):
		print('input hold_tree_file not found.')
		sys.exit(1)

	workdir = os.path.abspath(output_file_path)
	if not os.path.exists(workdir):
		os.mkdir(workdir)
	log = f"{workdir}/log"
	if not os.path.exists(log):
		os.mkdir(log)
	tmp = f"{workdir}/tmp"
	if not os.path.exists(tmp): 
		os.mkdir(tmp)
	err_model = f"{workdir}/err_model"
	if not os.path.exists(err_model):
		os.mkdir(err_model)
	k_o = f"{workdir}/k_output"
	if not os.path.exists(k_o):
		os.mkdir(k_o)
	k_p_log = f"{log}/k_output"
	if not os.path.exists(k_p_log):
		os.mkdir(k_p_log)
	k_t = f"{tmp}/k_output"
	if not os.path.exists(k_t):
		os.mkdir(k_t)


	cafe5_auto(od,pd,td,output_file_path,sp_name,hold_tree_file)



#!/nfs_genome/anaconda3/envs/metawrap-env/bin/python
# -*- coding: utf-8 -*-
from __future__ import print_function
####
## 1. 作者:	Shawn Yan shawn.yan@msn.com +86-15953380032
## 2. 版本:	0.03
## 3. 历史:	v0.01@2020/01/05; v0.02@20200114; v0.03@20200304
## 4. 目的:	根据conf和“sse05_bin.8_keggchart.csv”文件的定义，在conf对应的png文件上画矩形和连线，然后保存到“原名_color.png”中
## 5. 用法:	需要指定两个参数：1/2：conf文件的路径；2/2：“sse05_bin.8_keggchart.csv”文件的路径
## 注意：png文件的路径： 由conf计算出来
## 例如：	python kegg_chart.py conf文件路径 sse05_bin.8_keggchart.csv
## 5.1 单个文件执行:		python kegg_chart.py "归档2/map00920.conf"	sse05_bin.8_keggchart.csv
## 5.2 cmd命令行:		for %i in (*.conf) do python kegg_chart.py %i sse05_bin.8_keggchart.csv >>_run.log
## 5.3 将conf和png放在maps目录中： for %i in (maps/*.conf) do python kegg_chart.py %i sse05_bin.8_keggchart.csv >>_run.log
## 5.4 批处理run.bat:	for %%i in (*.conf) do python kegg_chart.py %%i sse05_bin.8_keggchart.csv >>_run.log
# 比如，run.bat，将不同运行批次的记录写到相应的log文件中，测试成功：
#set filename=%date:~0,4%%date:~5,2%%date:~8,2%_%time:~0,2%%time:~3,2%%time:~6,2%.log
#set "filename=%filename: =0%"
#rem echo %filename%
#for %%i in (maps\*.conf) do python kegg_chart_202003.py %%i sse05_bin.8_keggchart.csv >>_%filename%
####

import re			# 导入正则表达式模块
import cv2 as cv		#导入cv模块
import numpy as np	# 导入numpy模块
import sys			# 导入sys模块，用以接收python运行参数
import os

draw_rect_array = []	# 读取conf文件，再根据csv筛选出来的要进行绘制的列表，[0]为ko值，[1]为四个整数代表两个顶点
draw_line_array = []	# conf文件中的以line为开头的行中，若找到匹配的ko，则画线。第三版中增加的功能
global CODE_PATH	# 存储当前处理的PATHWAY通路名，即mapXXXXX.conf的koxxxxx
global CODE_KO	# 存储当前处理的ko值，即kXXXXX，实际没用到
global g_previous_kegg_key	# 用于排序后的draw_rect_array中间判断连续的两个项的ko是否相同，相同则画线
global g_previous_rect_set	# 同上

global matched_count

g_previous_kegg_key = ''			# 初始化
g_previous_rect_set = [0, 0, 0, 0]	# 初始化
matched_count = 0
# 4种颜色，对应0，1，2，3的4个val值
global CODE_COLOR
CODE_COLOR=[ # B G R
		[255, 127, 127],	# 蓝色1->浅蓝
		[000, 255, 255],	# 黄色2
		[000, 127, 255],	# 橙色3
		[000, 000, 255]	# 红色大于3的，没有就是白色
	]

# 从CSV中读取的二维字典
global PATHWAY_KO_COUNT
PATHWAY_KO_COUNT={}

# 二维Python字典，即散列表
# 本程序中，由第二个运行参数，即“sse05_bin.8_keggchart.csv”读取而来
# 第一列是key_a，第二列是key_b，第三列规整到<=3，为val
def ADD_TWO_DIM_DICT(thedict, key_a, key_b, val): 
	if key_a in thedict:
		thedict[key_a].update({key_b: val})
	else:
		thedict.update({key_a:{key_b: val}})



def DRAW_RECT (mat_img, kegg, rect):
	global g_previous_kegg_key
	global g_previous_rect_set

	# 1. 输出基本信息
	print("DRAW_RECT: PATHWAY_KO_COUNT[%s][%s]==%s" % (CODE_PATH, kegg, PATHWAY_KO_COUNT[CODE_PATH][kegg]))

	# 2. 与上一个的kegg值相同，画根线，好辩认？
	if kegg == g_previous_kegg_key: 
		cv.line(mat_img, (g_previous_rect_set[0]-1, g_previous_rect_set[1]-1), (rect[0]-1, rect[1]-1), 
				CODE_COLOR[int(PATHWAY_KO_COUNT[CODE_PATH][kegg])-1], 1, cv.LINE_AA)
		print ("Draw line for same kegg_key:", kegg, '	', end='')

	# 3. 先将区域内的白色像素，画上颜色。
	for xx in range(rect[0]+2, rect[2]+1):
		for yy in range(rect[1], rect[3]):
			if (not (mat_img[yy, xx] - [255, 255, 255]).all()):
				mat_img[yy, xx] = CODE_COLOR[int(PATHWAY_KO_COUNT[CODE_PATH][kegg])-1]

	# 4. 将第一列画上紫色，出头2像素，最后画这个线，避免被掩盖
	for yy in range(rect[1]-3, rect[3]+3):
		if (not (mat_img[yy, rect[0]+1] - [255, 255, 255]).all()):
			mat_img[yy, rect[0]+1] = [218, 112, 214] #淡紫色

	# 5. 更新全局变量
	g_previous_kegg_key = kegg
	g_previous_rect_set = rect
# END of def DRAW_RECT

## **** 函数定义结束 **** ##





if __name__ == '__main__':
	
	# 00. 读取程序的参数
	# 参数样例:CONF_File_Name='maps/map00920.conf' IMAGE_File_Name='maps/map00920.png'
	CONF_File_Name  = sys.argv[1]
	IMAGE_File_Name = CONF_File_Name.replace(".conf", ".png")
	print("CONF_File_Name == ", CONF_File_Name, ";\tIMAGE_File_Name == ", IMAGE_File_Name)
	
	CODE_PATH = 'M' + CONF_File_Name[len(CONF_File_Name)-10 : len(CONF_File_Name)-5]
	print("CODE_PATH == ", CODE_PATH)
	
	csvPath = sys.argv[2]


	# 01. 先读取并遍历CSV文件	#csvPath=="sse05_bin.8_keggchart.csv", 
	npCSV = np.loadtxt(csvPath, dtype = np.str, delimiter = "	")
	# csv文件每行遍历，将每一行的数据加载到“PATHWAY_KO_COUNT”二维数组中：
	for i in range(0, len(npCSV)):	
                #print(f"CODE_PATH:{CODE_PATH} {npCSV[i][0]}")
		if npCSV[i][0] == CODE_PATH:
			# 将“sse05_bin.8_keggchart.csv”文件中的第三列，强制转换为int整数类型
			# 再查查此字段的含义是：COUNT
			npCSV[i][2] = int(npCSV[i][2])
			if int(npCSV[i][2]) > 3:
				npCSV[i][2] = 3
                        
			ADD_TWO_DIM_DICT(PATHWAY_KO_COUNT, npCSV[i][0], npCSV[i][1], npCSV[i][2])
			#print(PATHWAY_KO_COUNT, npCSV[i][0], npCSV[i][1], npCSV[i][2])
	
	print("PATHWAY_KO_COUNT == ", PATHWAY_KO_COUNT)
	#### csv文件遍历结束 ####



	# 02. 打开conf文件，并遍历
	draw_rect_array = [] # 重新置零，要是以后在本程序批量处理的话，能避坑。
	draw_line_array = [] # 画线

	CONF_FILE_IN = open(CONF_File_Name, 'r')	# 传统行读取方式打开conf文件
	
	line_index = 0
	rect_found = 0

	# 承接约230行：print("%s[%s]=%s; "...
	print("PATHWAY_KO_COUNT matched: ")#, end=''
	try:
		while True:
			text_line = CONF_FILE_IN.readline()	# 读取conf文件的一行
			line_index += 1	# 行计数加1

			if not text_line: # 若读取行失败，则跳至下一个循环
				print("\nEND of CONF file read lines.", text_line)
				break	# 不要用continue，无法退出while循环

			# 行读取成功，则：
			#print("Line", line_index, type(text_line), text_line) # 输出行类型和行内容
			#if text_line[0:4] == 'rect': #OK 1 # 判断方法1
			
			# 若行不是以rect或line开头的，则跳至下一个while循环
			if not text_line.startswith('rect') and not text_line.startswith('line'):
				#print("line %d not started with rect or line: %s" % (line_index, text_line))
				continue
			 #OK 2# 判断方案2
			#print("Start with rect: ", text_line, end='')#检查当前行


				
			text_line = text_line.strip()	# 去掉末尾的 \r\n
			
			##---- 符合条件，正式执行的代码 ---- ##
			# 001: 当前行以制表符分隔，实际仅需[0] [1];
			# 样例：rect (281,244) (327,261)制表符/dbget-bin/www_bget?K04091+K00299+1.14.14.5+R07210
			# 接上行：制表符K04091 (ssuD), K00299 (ssuE), 1.14.14.5, R07210
			line_array = text_line.split("\t")
			#print(line_array[0] + " ---- " + line_array[1])# 输出检查行第一层分割结果

			# 002: RECT 再分割矩形定义字符串, 从第五个字符开始。OK
			line_array[0] = line_array[0].replace(') 1', '') #以“line (”开头的行，\t前有") 1"需要去掉
			line_array[0] = line_array[0][5:].replace('(', '').replace(')', '').replace(',', ' ')
			rect_array = line_array[0].split(' ')
			# 强制转化为整数类型
			for dd in range(0, len(rect_array)):
				rect_array[dd] = int(rect_array[dd])
				#rect_array[1] = int(rect_array[1])
				#rect_array[2] = int(rect_array[2])
				#rect_array[3] = int(rect_array[3])
			# 在次也就不要再截断了，有多少要多少，通通全收，适应line类型的多个数值
			#rect_array = rect_array[0:4] 
			#print("\n\nConf line Type: %s; rect_array: " % text_line[0:4], rect_array) #OK


			# 003: 再使用正则表达式分割第二列，例如：
			# /dbget-bin/www_bget?K00651+K00641+2.3.1.46+R01777
			matchObj = re.match( r'^(.*?\?)(.*)$', line_array[1], re.M|re.I)
			if matchObj:
# 001) 先获取问号后面的部分，即：matchObj.group(2)
#matchObj.group()==line_array[1];matchObj.group(1)=/dbget-bin/www_bget?
#matchObj.group(2)=K00651+K00641+2.3.1.46+R01777
				line_array[1] = matchObj.group(2)

# 002) 使用加号分割字符串
				ko_array_this_line = line_array[1].split('+')
				#print("len(ko_array_this_line)", len(ko_array_this_line), "ko_array_this_line", ko_array_this_line)
				
# 003) 分割后的各项，若不是以K或R开头的，则删除。
				for i in range(len(ko_array_this_line)-1, 0, -1):
					#print ("%d ko=%s" %(i, ko_array_this_line[i]))
					if not ko_array_this_line[i][0] in ['K', 'R']:
						del ko_array_this_line[i]

# 004) 再次类似第三项的检查：
#	删除一行仅一个类似KO值的情况，即未split成数组的情况
				if not ko_array_this_line[0][0] in ['K', 'R']:
					del ko_array_this_line[0]

# 005) 找完了第二列中的K和R，存储至ko_array_this_line，此时尚未区分rect和line两种类型的行
				#print("line 199 ko_array_this_line: ", ko_array_this_line)#, end = ''


# 此行中找到的K和R，再在csv中的定义中进行检索，
				ko_array_this_line_draw = []
				draw_line_by_this_line = False
# 找出CSV文件中定义的KO，没定义的略过


				for i in range(0, len(ko_array_this_line)):
# CODE_PATH为ko再加上5位数字，例如：CODE_PATH== ko01230
# 分别进行1级和2级检查：
# 	1. 此CODE_PATH在“sse05_bin.8_keggchart.csv”中是否有定义。
#		在第三个版本中，仅将csv中此ko值的数据填入PATHWAY_KO_COUNT中。
# 	2. 1成立的话，看此定义中是否有此ko_array_this_line[i]的定义。
# 1和2都成立的话，加入待画矩形的列表中。
					if CODE_PATH in PATHWAY_KO_COUNT \
						and ko_array_this_line[i] in PATHWAY_KO_COUNT[CODE_PATH]:
						if text_line.startswith('rect'):	# 若是rect行，则将此ko加入到此个矩形的绘制当中
							ko_array_this_line_draw.append(ko_array_this_line[i])
						else:	# 若是line开头的行，则进行标记，无需分拆矩形
							draw_line_by_this_line = True
							print("%s[%s]=%s; "%(text_line[0:4].upper(), ko_array_this_line[i], PATHWAY_KO_COUNT[CODE_PATH][ko_array_this_line[i]]), end='')
							# 在此输出csv中找的定义，以此判断是否需要画线
							#print("PATHWAY_KO_COUNT[%s]==%s" \
							#	% (ko_array_this_line[i], PATHWAY_KO_COUNT[CODE_PATH][ko_array_this_line[i]]))
						matched_count+=1
						if matched_count%5 == 0:
							print("")
				
				
				#print ("ko_array_this_line[%d] -> [%d]. " %(len(ko_array_this_line), len(ko_array_this_line_draw)))
				#print ("ko_array_this_line==", ko_array_this_line)
				if draw_line_by_this_line:
					draw_line_array.append(rect_array)
					#print("draw_line_array.append(rect_array)==", rect_array)
					#print("Line %d has line to draw: %s" %(line_index, text_line))
				
				# 再反向赋值，后面的代码尽量少改。
				ko_array_this_line = ko_array_this_line_draw
				if len(ko_array_this_line) > 0:
					#print(rect_array, "==", ko_array_this_line)
					rect_found += 1
					
					# v2 将一整个矩形划分成多个小矩形的单个宽度；累计后肯定会有结尾误差。
					perX = int((rect_array[2]-rect_array[0])/len(ko_array_this_line))
					# perY = int((rect_array[3]-rect_array[1])/len(ko_array_this_line)) #Y不需要分割
					for i in range(0, len(ko_array_this_line)): # 不要减一
						# 将conf中定义的矩形分割成条形码，每一条的第一列标紫。
						draw_rect_array.append([ko_array_this_line[i], 
							[rect_array[0] + i * perX,		rect_array[1], 
							 rect_array[0] + (i + 1) * perX,	rect_array[3]]
							])
				#else:	# end of: if len(ko_array_this_line) > 0:
					#print("Line %d has no rect to draw: %s" % (line_index, text_line))#, end=''

			else:
				print("Col 2 no match. ")#, end=''
			# end of: if matchObj:
		# end of: while True:
		print("PATHWAY_KO_COUNT matched end.")
	finally:
		CONF_FILE_IN.close()

	#print("draw_rect_array ==", draw_rect_array) #OK

	# 03. 开始绘图。
	draw_rect_array.sort() # 必须先排序，才能方便的在相同kegg值间画线！
	print("Rects to draw: %d;  Lines to draw: %d" % (len(draw_rect_array), len(draw_line_array) ))
	#print("draw_line_array == ", draw_line_array)
	if len(draw_rect_array) > 0 or len(draw_line_array) > 0:
		# 采用np读取图像，支持bmp、jpg、png、tiff等常用格式. 因为OpenCV无法读取中文名的图！
		# 例如：mat_img = cv.imread('归档2/map00920.png', cv.IMREAD_COLOR) # Failed due to Chinese charactors.
		mat_img = cv.imdecode(np.fromfile(IMAGE_File_Name, dtype=np.uint8), cv.IMREAD_COLOR)
		if mat_img is None:
			print("!!!!ERROR!!!! Image Read Failed: ", IMAGE_File_Name)
		else:
			print("Image Read Succeed: ", IMAGE_File_Name)
			#print("draw_rect_array", draw_rect_array)

			# 多个KO会使用同一个rect矩形定义，已经将矩形拆分成条形码方式。
			for da in range(0, len(draw_rect_array)):
				DRAW_RECT(mat_img, draw_rect_array[da][0], draw_rect_array[da][1])
				# 注：draw_rect_array定义：[[ko_code, [x1,y1,x2,y2]], ... ]
				
			font = cv.FONT_HERSHEY_SIMPLEX
			for dl in range(0, len(draw_line_array)):
				cv.putText(mat_img, 
					str(dl), # 要显示的文字
					(draw_line_array[dl][0] + 2, draw_line_array[dl][1] - 2), # 坐标，右上移动2像素
					font,		# 字体对象
					0.6,			# 字号放大
					(0, 0, 255),	# 颜色
					1			# 线条粗细
				)
				# 先统一向右上偏移2个像素
				for dt in range(0, len(draw_line_array[dl])):
					if dt %2==0:
						# X值向右偏移几个像素，例如2个，向右为+操作；
						draw_line_array[dl][dt]+=1
					else:
						# Y值向上偏移几个像素，例如2个，向上为-操作；
						draw_line_array[dl][dt]-=1
				for dt in range(0, int(len(draw_line_array[dl])/2)-1):
					cv.line(mat_img, 
						(draw_line_array[dl][dt*2+0], draw_line_array[dl][dt*2+1]), 
						(draw_line_array[dl][dt*2+2], draw_line_array[dl][dt*2+3]), 
						[255,0,255],
						#[(255-dl*4)%255,(0)%255, (255-dl*8)%255], #[255,0,255]这个颜色是紫色。
						1, # 线条粗细
						cv.LINE_8)	#cv.LINE_AA==反锯齿；在此反应过程一般是直线+拐角处的一些圆角，不要用LINE_AA，否则放大后比较难看。
			location = os.getcwd()	
			#save_path = IMAGE_File_Name.replace('.png', '_color.png')
			save_path = ''.join([location,"/",re.findall(r".*\/(M.*png)",IMAGE_File_Name)[0]])
			cv.imencode('.png', mat_img)[1].tofile(save_path)
			print("\nImage write done: ", save_path)
	else:
		print("\n!!!! Empty draw_rect_array & draw_line_array !!!!\nNo image file was not opened and saved:", IMAGE_File_Name)


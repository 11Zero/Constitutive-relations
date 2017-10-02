- abaqus中需要定义混凝土损伤破坏模型，因此查规范作出一系列本构关系曲线，基于GB2010混规
- matlab r2011b下的代码
- 其中C60-C80的这部分的混凝土变异系数在规范中未给出经验值，本人又没有试验数据，就取相邻强度近似值14%
- 对于混凝土抗拉抗压强度代表值，规范规定视情况决定取用标准值或设计值或平均值，这里都取了设计值
- 代码中：
- epsilon\_length代表模型选区的应变范围，为该数值乘以峰值应变
- epsilon\_step代表模型选区取值点每单位应变的密度，密度为每个峰值应变区间被划分的微段个数
- 导出xls文件
- 通过vbs脚本合并两个xls文件
- 通过vbs脚本对每个等级混凝土做应力应变曲线
- 最终生成concrete.xls

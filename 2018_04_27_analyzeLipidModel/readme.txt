20180427
ここでやっていることは，pptファイル
2018_04_24_報告用資料
2018_04_20_報告用資料
を見ると概略がわかる．
使用したcodeは，
"\\Azlab-fs01\東研究室\個人work\竹内(ひ)\codes\make_medium\forUSCTSim\makeRandomMedium.m"を参照のこと．
IMCL占有率の決定方法が明らかに間違っていることが判明している．
魚拓としてmakeRandomMedium_false.mを残しておく．
詳しくは，"\\Azlab-fs01\東研究室\個人work\竹内(ひ)\data\kwave\script\2018_09_28_realisticScatter_variousIMCL\readme.txt"を参照．


詳細
　グリッドサイズ
　　0.125 mm
　グリッド数
　　400うちセンサが占めるグリッド数320
　音源数
　　100(total: 200)
　関数作成
　　makeRandomMedium
　　getAverageSoundSpeed
　データセット
　　EMCLの占める割合
　　　1,2,3,4,5%
　　EMCLの個数
　　　1,2,3,4,5個
　　IMCLの割合
　　　0.2,0.4,0.6,0.8,1.0%
　　データセット数
　　　{5,3,0.2}を中心として直交表を作成．
　　　5+4+4の13通り
　AIC
　　halfWindowSize
　　　75
　　　　2018_03_05_中間試問スライドの42ページ参照
　　　　経験的に良さそうというだけ．
　  get_tof_AIC_from_singleRFData
　　　このexpのために作成．AveSSを出すのに有効．また，TTM (notTTDM)作成に有効．

　データセット追加
　　EMCLの占める割合
　　　6,8,10,12,14,16,18,20%と，大きく振ることで傾向を見てみる．
　　　（2018-04-30）

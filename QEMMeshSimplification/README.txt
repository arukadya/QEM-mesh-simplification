OS 		macOS Sonoma 14.4

コンパイラ	Apple clang version 15.0.0

コンパイルコマンド	g++ main.cpp Mesh.cpp -I <eigen3.4.0のパス> -std=c++17 -O2 -Wall

実行コマンド	./a.out

プログラムの概要	辺縮約によるメッシュ簡略化を行う．
		カレントディレクトリにあるlucy.objを読み込み，簡略化後の頂点リストと面リストをQEM_simplificate_lucy.objに出力する．

		Linear Solveは最小ノルム型(特異値分解とMoore-Penrose一般逆行列)を用いた．
		QEMの重みつけは面積重みで行った．
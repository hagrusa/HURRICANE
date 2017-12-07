make
python fabricate_data.py fabricate 500
./nbodymain data.txt .001 500 .1 4 RK4 0.2 result/
make clean
python fabricate_data.py
rm result/*.txt
ffmpeg -i result/graphic_result_%03d.png -pix_fmt yuv420p -r 8 simulation_viz.mp4
rm result/*.png

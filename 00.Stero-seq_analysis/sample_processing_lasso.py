import warnings
warnings.filterwarnings('ignore')
import stereo as st

input_gef_file = /path/to/tissue_gef_file
output_gem_file=/path/to/save_gem_file

input_image = /path/to/registered_tif
output_image = /path/to/selected_tif 

# read GEF file in BGEF
data = st.io.read_gef(input_gef_file, bin_type='bins', bin_size=10)

data.tl.cal_qc()ins = data.plt.interact_spatial_scatter(width=500, height=500, poly_select=True)
ins.show()

ins.export_high_res_area(
    input_gef_file,
    output_gem_file,
    # drop=True
)

ins.export_roi_image(
    origin_file_path=input_iamge,
    output_path=output_image,
    # drop=True
)


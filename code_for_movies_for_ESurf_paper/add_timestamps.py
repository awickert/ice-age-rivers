import glob
from PIL import Image
from PIL import ImageFont
from PIL import ImageDraw

model = 'G12'
# 26 from L and 04.1 for all but ICE-3G and ANU

filenames = sorted(glob.glob('map_'+model+'*'))

for filename in filenames:
  age = filename[-10:-4]
  print age
  age_numeric = int(age)
  agestr = '%04.1f' %(age_numeric/1000.)+' ka'
  img = Image.open(filename)
  draw = ImageDraw.Draw(img)
  # font = ImageFont.truetype(<font-file>, <font-size>)
  #fontpath = '/usr/share/fonts/truetype/adf/GilliusADF-Bold.otf'
  #fontpath = '/usr/share/fonts/truetype/ubuntu-font-family/UbuntuMono-B.ttf'
  #fontpath = '/usr/share/fonts/truetype/droid/DroidSansMono.ttf'
  fontpath = '/usr/share/fonts/truetype/ubuntu-font-family/Ubuntu-R.ttf'
  font = ImageFont.truetype(fontpath, 60)
  # draw.text((x, y),"Sample Text",(r,g,b))
  top_bottom_offset = draw.textsize(agestr, font=font)[-1]
  draw.text((16, img.getbbox()[-1] - 26 - top_bottom_offset), \
             agestr,(255,255,255), font=font)
  img.save('timestamped/'+filename)

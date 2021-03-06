# Recorded script from Mayavi2
from numpy import array
try:
    engine = mayavi.engine
except NameError:
    from enthought.mayavi.api import Engine
    engine = Engine()
    engine.start()
if len(engine.scenes) == 0:
    engine.new_scene()
# ------------------------------------------- 
scene = engine.scenes[0]
scene.scene.background = (0.49803921568627452, 0.49803921568627452, 0.49803921568627452)
scene.scene.background = (0.97254901960784312, 0.97254901960784312, 0.97254901960784312)
module_manager = engine.scenes[0].children[0].children[0]
module_manager.scalar_lut_manager.title_text_property.shadow_offset = array([ 1, -1])
module_manager.scalar_lut_manager.title_text_property.color = (0.039215686274509803, 0.039215686274509803, 0.039215686274509803)
module_manager.scalar_lut_manager.label_text_property.shadow_offset = array([ 1, -1])
module_manager.scalar_lut_manager.label_text_property.color = (0.039215686274509803, 0.039215686274509803, 0.039215686274509803)
module_manager.vector_lut_manager.title_text_property.shadow_offset = array([ 1, -1])
module_manager.vector_lut_manager.title_text_property.color = (0.039215686274509803, 0.039215686274509803, 0.039215686274509803)
module_manager.vector_lut_manager.label_text_property.shadow_offset = array([ 1, -1])
module_manager.vector_lut_manager.label_text_property.color = (0.039215686274509803, 0.039215686274509803, 0.039215686274509803)
vectors = engine.scenes[0].children[0].children[0].children[0]
vectors.actor.property.specular_color = (0.039215686274509803, 0.039215686274509803, 0.039215686274509803)
vectors.actor.property.diffuse_color = (0.039215686274509803, 0.039215686274509803, 0.039215686274509803)
vectors.actor.property.ambient_color = (0.039215686274509803, 0.039215686274509803, 0.039215686274509803)
vectors.actor.property.color = (0.039215686274509803, 0.039215686274509803, 0.039215686274509803)
text = engine.scenes[0].children[0].children[0].children[1]
text.property.shadow_offset = array([ 1, -1])
text.property.color = (0.039215686274509803, 0.039215686274509803, 0.039215686274509803)
text1 = engine.scenes[0].children[0].children[0].children[2]
text1.property.shadow_offset = array([ 1, -1])
text1.property.color = (0.039215686274509803, 0.039215686274509803, 0.039215686274509803)
text2 = engine.scenes[0].children[0].children[0].children[3]
text2.property.shadow_offset = array([ 1, -1])
text2.property.color = (0.039215686274509803, 0.039215686274509803, 0.039215686274509803)
text3 = engine.scenes[0].children[0].children[0].children[4]
text3.property.shadow_offset = array([ 1, -1])
text3.property.color = (0.039215686274509803, 0.039215686274509803, 0.039215686274509803)
text4 = engine.scenes[0].children[0].children[0].children[5]
text4.property.shadow_offset = array([ 1, -1])
text4.property.color = (0.039215686274509803, 0.039215686274509803, 0.039215686274509803)
text5 = engine.scenes[0].children[0].children[0].children[6]
text5.property.shadow_offset = array([ 1, -1])
text5.property.color = (0.039215686274509803, 0.039215686274509803, 0.039215686274509803)
text6 = engine.scenes[0].children[0].children[0].children[7]
text6.property.shadow_offset = array([ 1, -1])
text6.property.color = (0.039215686274509803, 0.039215686274509803, 0.039215686274509803)
text7 = engine.scenes[0].children[0].children[0].children[8]
text7.property.shadow_offset = array([ 1, -1])
text7.property.color = (0.039215686274509803, 0.039215686274509803, 0.039215686274509803)
text8 = engine.scenes[0].children[0].children[0].children[9]
text8.property.shadow_offset = array([ 1, -1])
text8.property.color = (0.039215686274509803, 0.039215686274509803, 0.039215686274509803)
text9 = engine.scenes[0].children[0].children[0].children[10]
text9.property.shadow_offset = array([ 1, -1])
text9.property.color = (0.039215686274509803, 0.039215686274509803, 0.039215686274509803)
text10 = engine.scenes[0].children[0].children[0].children[11]
text10.property.shadow_offset = array([ 1, -1])
text10.property.color = (0.039215686274509803, 0.039215686274509803, 0.039215686274509803)
text11 = engine.scenes[0].children[0].children[0].children[12]
text11.property.shadow_offset = array([ 1, -1])
text11.property.color = (0.039215686274509803, 0.039215686274509803, 0.039215686274509803)
text12 = engine.scenes[0].children[0].children[0].children[13]
text12.property.shadow_offset = array([ 1, -1])
text12.property.color = (0.039215686274509803, 0.039215686274509803, 0.039215686274509803)
text13 = engine.scenes[0].children[0].children[0].children[14]
text13.property.shadow_offset = array([ 1, -1])
text13.property.color = (0.039215686274509803, 0.039215686274509803, 0.039215686274509803)
text14 = engine.scenes[0].children[0].children[0].children[15]
text14.property.shadow_offset = array([ 1, -1])
text14.property.color = (0.039215686274509803, 0.039215686274509803, 0.039215686274509803)
text15 = engine.scenes[0].children[0].children[0].children[16]
text15.property.shadow_offset = array([ 1, -1])
text15.property.color = (0.039215686274509803, 0.039215686274509803, 0.039215686274509803)
text16 = engine.scenes[0].children[0].children[0].children[17]
text16.property.shadow_offset = array([ 1, -1])
text16.property.color = (0.039215686274509803, 0.039215686274509803, 0.039215686274509803)
text17 = engine.scenes[0].children[0].children[0].children[18]
text17.property.shadow_offset = array([ 1, -1])
text17.property.color = (0.039215686274509803, 0.039215686274509803, 0.039215686274509803)
text18 = engine.scenes[0].children[0].children[0].children[19]
text18.property.shadow_offset = array([ 1, -1])
text18.property.color = (0.039215686274509803, 0.039215686274509803, 0.039215686274509803)
text19 = engine.scenes[0].children[0].children[0].children[20]
text19.property.shadow_offset = array([ 1, -1])
text19.property.color = (0.039215686274509803, 0.039215686274509803, 0.039215686274509803)
text21 = engine.scenes[0].children[0].children[0].children[22]
text21.property.shadow_offset = array([ 1, -1])
text21.property.color = (0.039215686274509803, 0.039215686274509803, 0.039215686274509803)
text22 = engine.scenes[0].children[0].children[0].children[23]
text22.property.shadow_offset = array([ 1, -1])
text22.property.color = (0.039215686274509803, 0.039215686274509803, 0.039215686274509803)
module_manager1 = engine.scenes[0].children[1].children[0].children[0].children[0]
module_manager1.scalar_lut_manager.title_text_property.shadow_offset = array([ 1, -1])
module_manager1.scalar_lut_manager.title_text_property.color = (0.039215686274509803, 0.039215686274509803, 0.039215686274509803)
module_manager1.scalar_lut_manager.label_text_property.shadow_offset = array([ 1, -1])
module_manager1.scalar_lut_manager.label_text_property.color = (0.039215686274509803, 0.039215686274509803, 0.039215686274509803)
module_manager1.vector_lut_manager.title_text_property.shadow_offset = array([ 1, -1])
module_manager1.vector_lut_manager.title_text_property.color = (0.039215686274509803, 0.039215686274509803, 0.039215686274509803)
module_manager1.vector_lut_manager.label_text_property.shadow_offset = array([ 1, -1])
module_manager1.vector_lut_manager.label_text_property.color = (0.039215686274509803, 0.039215686274509803, 0.039215686274509803)
surface = engine.scenes[0].children[1].children[0].children[0].children[0].children[0]
surface.actor.property.specular_color = (0.039215686274509803, 0.039215686274509803, 0.039215686274509803)
surface.actor.property.diffuse_color = (0.039215686274509803, 0.039215686274509803, 0.039215686274509803)
surface.actor.property.ambient_color = (0.039215686274509803, 0.039215686274509803, 0.039215686274509803)
surface.actor.property.color = (0.039215686274509803, 0.039215686274509803, 0.039215686274509803)
scene.scene.foreground = (0.039215686274509803, 0.039215686274509803, 0.039215686274509803)
module_manager.scalar_lut_manager.title_text_property.shadow_offset = array([ 1, -1])
module_manager.scalar_lut_manager.title_text_property.color = (0.027450980392156862, 0.027450980392156862, 0.027450980392156862)
module_manager.scalar_lut_manager.label_text_property.shadow_offset = array([ 1, -1])
module_manager.scalar_lut_manager.label_text_property.color = (0.027450980392156862, 0.027450980392156862, 0.027450980392156862)
module_manager.vector_lut_manager.title_text_property.shadow_offset = array([ 1, -1])
module_manager.vector_lut_manager.title_text_property.color = (0.027450980392156862, 0.027450980392156862, 0.027450980392156862)
module_manager.vector_lut_manager.label_text_property.shadow_offset = array([ 1, -1])
module_manager.vector_lut_manager.label_text_property.color = (0.027450980392156862, 0.027450980392156862, 0.027450980392156862)
vectors.actor.property.specular_color = (0.027450980392156862, 0.027450980392156862, 0.027450980392156862)
vectors.actor.property.diffuse_color = (0.027450980392156862, 0.027450980392156862, 0.027450980392156862)
vectors.actor.property.ambient_color = (0.027450980392156862, 0.027450980392156862, 0.027450980392156862)
vectors.actor.property.color = (0.027450980392156862, 0.027450980392156862, 0.027450980392156862)
text.property.shadow_offset = array([ 1, -1])
text.property.color = (0.027450980392156862, 0.027450980392156862, 0.027450980392156862)
text1.property.shadow_offset = array([ 1, -1])
text1.property.color = (0.027450980392156862, 0.027450980392156862, 0.027450980392156862)
text2.property.shadow_offset = array([ 1, -1])
text2.property.color = (0.027450980392156862, 0.027450980392156862, 0.027450980392156862)
text3.property.shadow_offset = array([ 1, -1])
text3.property.color = (0.027450980392156862, 0.027450980392156862, 0.027450980392156862)
text4.property.shadow_offset = array([ 1, -1])
text4.property.color = (0.027450980392156862, 0.027450980392156862, 0.027450980392156862)
text5.property.shadow_offset = array([ 1, -1])
text5.property.color = (0.027450980392156862, 0.027450980392156862, 0.027450980392156862)
text6.property.shadow_offset = array([ 1, -1])
text6.property.color = (0.027450980392156862, 0.027450980392156862, 0.027450980392156862)
text7.property.shadow_offset = array([ 1, -1])
text7.property.color = (0.027450980392156862, 0.027450980392156862, 0.027450980392156862)
text8.property.shadow_offset = array([ 1, -1])
text8.property.color = (0.027450980392156862, 0.027450980392156862, 0.027450980392156862)
text9.property.shadow_offset = array([ 1, -1])
text9.property.color = (0.027450980392156862, 0.027450980392156862, 0.027450980392156862)
text10.property.shadow_offset = array([ 1, -1])
text10.property.color = (0.027450980392156862, 0.027450980392156862, 0.027450980392156862)
text11.property.shadow_offset = array([ 1, -1])
text11.property.color = (0.027450980392156862, 0.027450980392156862, 0.027450980392156862)
text12.property.shadow_offset = array([ 1, -1])
text12.property.color = (0.027450980392156862, 0.027450980392156862, 0.027450980392156862)
text13.property.shadow_offset = array([ 1, -1])
text13.property.color = (0.027450980392156862, 0.027450980392156862, 0.027450980392156862)
text14.property.shadow_offset = array([ 1, -1])
text14.property.color = (0.027450980392156862, 0.027450980392156862, 0.027450980392156862)
text15.property.shadow_offset = array([ 1, -1])
text15.property.color = (0.027450980392156862, 0.027450980392156862, 0.027450980392156862)
text16.property.shadow_offset = array([ 1, -1])
text16.property.color = (0.027450980392156862, 0.027450980392156862, 0.027450980392156862)
text17.property.shadow_offset = array([ 1, -1])
text17.property.color = (0.027450980392156862, 0.027450980392156862, 0.027450980392156862)
text18.property.shadow_offset = array([ 1, -1])
text18.property.color = (0.027450980392156862, 0.027450980392156862, 0.027450980392156862)
text19.property.shadow_offset = array([ 1, -1])
text19.property.color = (0.027450980392156862, 0.027450980392156862, 0.027450980392156862)
module_manager1.scalar_lut_manager.title_text_property.shadow_offset = array([ 1, -1])
module_manager1.scalar_lut_manager.title_text_property.color = (0.027450980392156862, 0.027450980392156862, 0.027450980392156862)
module_manager1.scalar_lut_manager.label_text_property.shadow_offset = array([ 1, -1])
module_manager1.scalar_lut_manager.label_text_property.color = (0.027450980392156862, 0.027450980392156862, 0.027450980392156862)
module_manager1.vector_lut_manager.title_text_property.shadow_offset = array([ 1, -1])
module_manager1.vector_lut_manager.title_text_property.color = (0.027450980392156862, 0.027450980392156862, 0.027450980392156862)
module_manager1.vector_lut_manager.label_text_property.shadow_offset = array([ 1, -1])
module_manager1.vector_lut_manager.label_text_property.color = (0.027450980392156862, 0.027450980392156862, 0.027450980392156862)
surface.actor.property.specular_color = (0.027450980392156862, 0.027450980392156862, 0.027450980392156862)
surface.actor.property.diffuse_color = (0.027450980392156862, 0.027450980392156862, 0.027450980392156862)
surface.actor.property.ambient_color = (0.027450980392156862, 0.027450980392156862, 0.027450980392156862)
surface.actor.property.color = (0.027450980392156862, 0.027450980392156862, 0.027450980392156862)
scene.scene.foreground = (0.027450980392156862, 0.027450980392156862, 0.027450980392156862)
module_manager.scalar_lut_manager.title_text_property.shadow_offset = array([ 1, -1])
module_manager.scalar_lut_manager.title_text_property.color = (0.011764705882352941, 0.011764705882352941, 0.011764705882352941)
module_manager.scalar_lut_manager.label_text_property.shadow_offset = array([ 1, -1])
module_manager.scalar_lut_manager.label_text_property.color = (0.011764705882352941, 0.011764705882352941, 0.011764705882352941)
module_manager.vector_lut_manager.title_text_property.shadow_offset = array([ 1, -1])
module_manager.vector_lut_manager.title_text_property.color = (0.011764705882352941, 0.011764705882352941, 0.011764705882352941)
module_manager.vector_lut_manager.label_text_property.shadow_offset = array([ 1, -1])
module_manager.vector_lut_manager.label_text_property.color = (0.011764705882352941, 0.011764705882352941, 0.011764705882352941)
vectors.actor.property.specular_color = (0.011764705882352941, 0.011764705882352941, 0.011764705882352941)
vectors.actor.property.diffuse_color = (0.011764705882352941, 0.011764705882352941, 0.011764705882352941)
vectors.actor.property.ambient_color = (0.011764705882352941, 0.011764705882352941, 0.011764705882352941)
vectors.actor.property.color = (0.011764705882352941, 0.011764705882352941, 0.011764705882352941)
text.property.shadow_offset = array([ 1, -1])
text.property.color = (0.011764705882352941, 0.011764705882352941, 0.011764705882352941)
text1.property.shadow_offset = array([ 1, -1])
text1.property.color = (0.011764705882352941, 0.011764705882352941, 0.011764705882352941)
text2.property.shadow_offset = array([ 1, -1])
text2.property.color = (0.011764705882352941, 0.011764705882352941, 0.011764705882352941)
text3.property.shadow_offset = array([ 1, -1])
text3.property.color = (0.011764705882352941, 0.011764705882352941, 0.011764705882352941)
text4.property.shadow_offset = array([ 1, -1])
text4.property.color = (0.011764705882352941, 0.011764705882352941, 0.011764705882352941)
text5.property.shadow_offset = array([ 1, -1])
text5.property.color = (0.011764705882352941, 0.011764705882352941, 0.011764705882352941)
text6.property.shadow_offset = array([ 1, -1])
text6.property.color = (0.011764705882352941, 0.011764705882352941, 0.011764705882352941)
text7.property.shadow_offset = array([ 1, -1])
text7.property.color = (0.011764705882352941, 0.011764705882352941, 0.011764705882352941)
text8.property.shadow_offset = array([ 1, -1])
text8.property.color = (0.011764705882352941, 0.011764705882352941, 0.011764705882352941)
text9.property.shadow_offset = array([ 1, -1])
text9.property.color = (0.011764705882352941, 0.011764705882352941, 0.011764705882352941)
text10.property.shadow_offset = array([ 1, -1])
text10.property.color = (0.011764705882352941, 0.011764705882352941, 0.011764705882352941)
text11.property.shadow_offset = array([ 1, -1])
text11.property.color = (0.011764705882352941, 0.011764705882352941, 0.011764705882352941)
text12.property.shadow_offset = array([ 1, -1])
text12.property.color = (0.011764705882352941, 0.011764705882352941, 0.011764705882352941)
text13.property.shadow_offset = array([ 1, -1])
text13.property.color = (0.011764705882352941, 0.011764705882352941, 0.011764705882352941)
text14.property.shadow_offset = array([ 1, -1])
text14.property.color = (0.011764705882352941, 0.011764705882352941, 0.011764705882352941)
text15.property.shadow_offset = array([ 1, -1])
text15.property.color = (0.011764705882352941, 0.011764705882352941, 0.011764705882352941)
text16.property.shadow_offset = array([ 1, -1])
text16.property.color = (0.011764705882352941, 0.011764705882352941, 0.011764705882352941)
text17.property.shadow_offset = array([ 1, -1])
text17.property.color = (0.011764705882352941, 0.011764705882352941, 0.011764705882352941)
text18.property.shadow_offset = array([ 1, -1])
text18.property.color = (0.011764705882352941, 0.011764705882352941, 0.011764705882352941)
text19.property.shadow_offset = array([ 1, -1])
text19.property.color = (0.011764705882352941, 0.011764705882352941, 0.011764705882352941)
module_manager1.scalar_lut_manager.title_text_property.shadow_offset = array([ 1, -1])
module_manager1.scalar_lut_manager.title_text_property.color = (0.011764705882352941, 0.011764705882352941, 0.011764705882352941)
module_manager1.scalar_lut_manager.label_text_property.shadow_offset = array([ 1, -1])
module_manager1.scalar_lut_manager.label_text_property.color = (0.011764705882352941, 0.011764705882352941, 0.011764705882352941)
module_manager1.vector_lut_manager.title_text_property.shadow_offset = array([ 1, -1])
module_manager1.vector_lut_manager.title_text_property.color = (0.011764705882352941, 0.011764705882352941, 0.011764705882352941)
module_manager1.vector_lut_manager.label_text_property.shadow_offset = array([ 1, -1])
module_manager1.vector_lut_manager.label_text_property.color = (0.011764705882352941, 0.011764705882352941, 0.011764705882352941)
surface.actor.property.specular_color = (0.011764705882352941, 0.011764705882352941, 0.011764705882352941)
surface.actor.property.diffuse_color = (0.011764705882352941, 0.011764705882352941, 0.011764705882352941)
surface.actor.property.ambient_color = (0.011764705882352941, 0.011764705882352941, 0.011764705882352941)
surface.actor.property.color = (0.011764705882352941, 0.011764705882352941, 0.011764705882352941)
scene.scene.foreground = (0.011764705882352941, 0.011764705882352941, 0.011764705882352941)
scene.scene.background = (0.9882352941176471, 0.9882352941176471, 0.9882352941176471)
scene.scene.background = (1.0, 1.0, 1.0)
module_manager.scalar_lut_manager.title_text_property.shadow_offset = array([ 1, -1])
module_manager.scalar_lut_manager.title_text_property.color = (0.0, 0.0, 0.0)
module_manager.scalar_lut_manager.label_text_property.shadow_offset = array([ 1, -1])
module_manager.scalar_lut_manager.label_text_property.color = (0.0, 0.0, 0.0)
module_manager.vector_lut_manager.title_text_property.shadow_offset = array([ 1, -1])
module_manager.vector_lut_manager.title_text_property.color = (0.0, 0.0, 0.0)
module_manager.vector_lut_manager.label_text_property.shadow_offset = array([ 1, -1])
module_manager.vector_lut_manager.label_text_property.color = (0.0, 0.0, 0.0)
vectors.actor.property.specular_color = (0.0, 0.0, 0.0)
vectors.actor.property.diffuse_color = (0.0, 0.0, 0.0)
vectors.actor.property.ambient_color = (0.0, 0.0, 0.0)
vectors.actor.property.color = (0.0, 0.0, 0.0)
text.property.shadow_offset = array([ 1, -1])
text.property.color = (0.0, 0.0, 0.0)
text1.property.shadow_offset = array([ 1, -1])
text1.property.color = (0.0, 0.0, 0.0)
text2.property.shadow_offset = array([ 1, -1])
text2.property.color = (0.0, 0.0, 0.0)
text3.property.shadow_offset = array([ 1, -1])
text3.property.color = (0.0, 0.0, 0.0)
text4.property.shadow_offset = array([ 1, -1])
text4.property.color = (0.0, 0.0, 0.0)
text5.property.shadow_offset = array([ 1, -1])
text5.property.color = (0.0, 0.0, 0.0)
text6.property.shadow_offset = array([ 1, -1])
text6.property.color = (0.0, 0.0, 0.0)
text7.property.shadow_offset = array([ 1, -1])
text7.property.color = (0.0, 0.0, 0.0)
text8.property.shadow_offset = array([ 1, -1])
text8.property.color = (0.0, 0.0, 0.0)
text9.property.shadow_offset = array([ 1, -1])
text9.property.color = (0.0, 0.0, 0.0)
text10.property.shadow_offset = array([ 1, -1])
text10.property.color = (0.0, 0.0, 0.0)
text11.property.shadow_offset = array([ 1, -1])
text11.property.color = (0.0, 0.0, 0.0)
text12.property.shadow_offset = array([ 1, -1])
text12.property.color = (0.0, 0.0, 0.0)
text13.property.shadow_offset = array([ 1, -1])
text13.property.color = (0.0, 0.0, 0.0)
text14.property.shadow_offset = array([ 1, -1])
text14.property.color = (0.0, 0.0, 0.0)
text15.property.shadow_offset = array([ 1, -1])
text15.property.color = (0.0, 0.0, 0.0)
text16.property.shadow_offset = array([ 1, -1])
text16.property.color = (0.0, 0.0, 0.0)
text17.property.shadow_offset = array([ 1, -1])
text17.property.color = (0.0, 0.0, 0.0)
text18.property.shadow_offset = array([ 1, -1])
text18.property.color = (0.0, 0.0, 0.0)
text19.property.shadow_offset = array([ 1, -1])
text19.property.color = (0.0, 0.0, 0.0)
module_manager1.scalar_lut_manager.title_text_property.shadow_offset = array([ 1, -1])
module_manager1.scalar_lut_manager.title_text_property.color = (0.0, 0.0, 0.0)
module_manager1.scalar_lut_manager.label_text_property.shadow_offset = array([ 1, -1])
module_manager1.scalar_lut_manager.label_text_property.color = (0.0, 0.0, 0.0)
module_manager1.vector_lut_manager.title_text_property.shadow_offset = array([ 1, -1])
module_manager1.vector_lut_manager.title_text_property.color = (0.0, 0.0, 0.0)
module_manager1.vector_lut_manager.label_text_property.shadow_offset = array([ 1, -1])
module_manager1.vector_lut_manager.label_text_property.color = (0.0, 0.0, 0.0)
surface.actor.property.specular_color = (0.0, 0.0, 0.0)
surface.actor.property.diffuse_color = (0.0, 0.0, 0.0)
surface.actor.property.ambient_color = (0.0, 0.0, 0.0)
surface.actor.property.color = (0.0, 0.0, 0.0)
scene.scene.foreground = (0.0, 0.0, 0.0)
module_manager.scalar_lut_manager.scalar_bar.position2 = array([ 0.8 ,  0.17])
module_manager.scalar_lut_manager.scalar_bar.reference_count = 4
module_manager.scalar_lut_manager.scalar_bar.position = array([ 0.1 ,  0.01])
module_manager.scalar_lut_manager.show_scalar_bar = True
module_manager.scalar_lut_manager.title_text_property.shadow_offset = array([ 1, -1])
module_manager.scalar_lut_manager.title_text_property.shadow = True
module_manager.scalar_lut_manager.label_text_property.shadow_offset = array([ 1, -1])
module_manager.scalar_lut_manager.label_text_property.shadow = True
module_manager.scalar_lut_manager.shadow = True
module_manager.scalar_lut_manager.use_default_name = False
module_manager.scalar_lut_manager.scalar_bar.position2 = array([ 0.8 ,  0.17])
module_manager.scalar_lut_manager.scalar_bar.position = array([ 0.1 ,  0.01])
module_manager.scalar_lut_manager.scalar_bar.title = 'Enrichment Ratio'
module_manager.scalar_lut_manager.scalar_bar.position2 = array([ 0.8 ,  0.17])
module_manager.scalar_lut_manager.scalar_bar.position = array([ 0.1 ,  0.01])
module_manager.scalar_lut_manager.data_name = u'Enrichment Ratio'
module_manager.scalar_lut_manager.use_default_range = False
module_manager.scalar_lut_manager.scalar_bar.position2 = array([ 0.8 ,  0.17])
module_manager.scalar_lut_manager.scalar_bar.position = array([ 0.1 ,  0.01])
module_manager.scalar_lut_manager.data_range = array([ 0.     ,  3.50446])
module_manager.scalar_lut_manager.scalar_bar.position2 = array([ 0.8 ,  0.17])
module_manager.scalar_lut_manager.scalar_bar.position = array([ 0.1 ,  0.01])
module_manager.scalar_lut_manager.data_range = array([ 0.     ,  3.50446])
module_manager.scalar_lut_manager.scalar_bar.position2 = array([ 0.8 ,  0.17])
module_manager.scalar_lut_manager.scalar_bar.position = array([ 0.1 ,  0.01])
module_manager.scalar_lut_manager.data_range = array([ 0. ,  3.3])
module_manager.scalar_lut_manager.scalar_bar.position2 = array([ 0.8 ,  0.17])
module_manager.scalar_lut_manager.scalar_bar.position = array([ 0.1 ,  0.01])
module_manager.scalar_lut_manager.data_range = array([ 0. ,  3.3])

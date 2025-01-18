import bpy
import addon_utils
import sys

addon_name = "batoms"

# Enable the addon
try:
    addon_utils.enable(addon_name, default_set=True)
    bpy.ops.wm.save_userpref()
    print(f"✅ Successfully enabled Blender addon: {addon_name}")
except Exception as e:
    print(f"❌ Failed to enable addon {addon_name}: {e}")
    sys.exit(1)

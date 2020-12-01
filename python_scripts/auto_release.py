# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#   Authors: T. Lemane
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU Affero General Public License as
#  published by the Free Software Foundation, either version 3 of the
#  License, or (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Affero General Public License for more details.
#
#  You should have received a copy of the GNU Affero General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

from github import Github
import os
import sys
import platform
import importlib

py_path = os.path.dirname(os.path.abspath(__file__))
sys.path.append(py_path+'/..')

system = platform.system()

LINUX, OSX = system == 'Linux', system == 'Darwin'

d = sys.argv[1]
package = ''.join([p for p in os.listdir(d) if p.endswith('.tar.gz')])
if not package:
    print('Package not found')
    sys.exit(0)

ori_version = ''.join(package.split('-')[1])
print(ori_version)
PACKAGE_NAME = f'ORI-{ori_version}-bin-{system}.tar.gz'
PACKAGE_PATH = f'{d}/{PACKAGE_NAME}'

GH_TOKEN = os.environ['GH_TOKEN']

BUILD_INFO = ''

g = Github(GH_TOKEN)
u = g.get_user()
ori_repo = u.get_repo('ORI')

LATEST_VERSION = []
try:
    latest_release = ori_repo.get_latest_release()
    assets = latest_release.get_assets()
    assets_name = [asset.name for asset in assets]
    tag_name = latest_release.tag_name
    tag_name = tag_name[1:] if tag_name.startswith('v') else tag_name
    LATEST_VERSION = list(map(int, tag_name.split('.')))
except:
    pass

CURRENT_VERSION = list(map(int, ori_version.split('.')))

print(BUILD_INFO)
print(system)
print(PACKAGE_PATH)
print(LATEST_VERSION, CURRENT_VERSION)

if not CURRENT_VERSION == [0,0,0]:
    if CURRENT_VERSION > LATEST_VERSION or not LATEST_VERSION:
        release = ori_repo.create_git_release(f'v{ori_version}', f'Release v{ori_version}', BUILD_INFO, False, prerelease=False)
        release.upload_asset(PACKAGE_PATH)
    elif CURRENT_VERSION == LATEST_VERSION:
        if PACKAGE_NAME not in assets_name:
            latest_release.upload_asset(PACKAGE_PATH)
    else:
        print('Nothing to upload')
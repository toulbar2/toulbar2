import os
import re
import sys
import platform
import subprocess

from setuptools import setup
from setuptools.extension import Extension
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
from distutils.version import LooseVersion

python_Path = sys.exec_prefix
print("python_Path", python_Path)

class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=''):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)


class CMakeBuild(build_ext):
    def run(self):
        try:
            out = subprocess.check_output(['cmake', '--version'])
        except OSError:
            raise RuntimeError("CMake must be installed to build the following extensions: " +
                               ", ".join(e.name for e in self.extensions))

        cmake_version = LooseVersion(re.search(r'version\s*([\d.]+)', out.decode()).group(1))
        if cmake_version < LooseVersion('3.5.0'):
            raise RuntimeError("CMake >= 3.5.0 is required")

        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))
        cmake_args = ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=' + extdir, '-DPython3_ROOT_DIR=' + python_Path,
                      '-DPYTB2=ON']

        build_type = os.environ.get("BUILD_TYPE", "Release")
        build_args = ['--config', build_type]

        # Pile all .so in one place and use $ORIGIN as RPATH
        cmake_args += ["-DCMAKE_BUILD_WITH_INSTALL_RPATH=TRUE"]
        cmake_args += ["-DCMAKE_INSTALL_RPATH={}".format("$ORIGIN")]

        if platform.system() == "Windows":
            cmake_args += ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{}={}'.format(build_type.upper(), extdir)]
            if sys.maxsize > 2**32:
                cmake_args += ['-A', 'x64']
            build_args += ['--', '/m']
        else:
            cmake_args += ['-DCMAKE_BUILD_TYPE=' + build_type]
            build_args += ['--', '-j4']

        env = os.environ.copy()
        env['CXXFLAGS'] = '{} -DVERSION_INFO=\\"{}\\"'.format(env.get('CXXFLAGS', ''),
                                                              self.distribution.get_version())
        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)
        subprocess.check_call(['cmake', ext.sourcedir] + cmake_args, cwd=self.build_temp, env=env)
        subprocess.check_call(['cmake',
                               '--build', '.',
                               '--target', os.path.basename('pytb2')
                               ] + build_args,
                              cwd=self.build_temp)

def read_version():
    return subprocess.run(['git','describe','--abbrev=0','--tags','--always'],capture_output=True).stdout.strip().decode("utf-8")

setup(
    name='pytoulbar2',
    version="0.0.0.2",
    author='ToulBar2 team',
    author_email='thomas.schiex@inrae.fr',
    description='ToulBar2 Python package',
    long_description_content_type="text/markdown",
    license='MIT',
    keywords='optimization graphical-model',
    long_description=open("README.md").read(),
    ext_modules=[CMakeExtension('pytoulbar2.pytb2')],
    packages=['pytoulbar2'],
    package_dir={'pytoulbar2': 'pytoulbar2'},
    cmdclass=dict(build_ext=CMakeBuild),
    url='http://miat.inrae.fr/toulbar2',
    zip_safe=False,
    classifiers=[
        'License :: OSI Approved :: MIT License',
        'Topic :: Scientific/Engineering :: Artificial Intelligence',
        'Topic :: Scientific/Engineering :: Mathematics',
        'Programming Language :: C++',
        'Programming Language :: Python :: 3',
        'Intended Audience :: Science/Research'
    ],
)


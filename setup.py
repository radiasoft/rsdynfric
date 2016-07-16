# -*- coding: utf-8 -*-
u"""rsfriction setup script

:copyright: Copyright (c) 2016 RadiaSoft LLC.  All Rights Reserved.
:license: http://www.apache.org/licenses/LICENSE-2.0.html
"""
try:
    from pykern import pksetup
except ImportError:
    import pip
    pip.main(['install', '--upgrade', 'pykern'])
    from pykern import pksetup


pksetup.setup(
    name='rsfriction',
    author='RadiaSoft LLC',
    author_email='pip@radiasoft.net',
    description='Python library for simulating the dynamical friction force exerted on ions by magnetized electrons',
    license='http://www.apache.org/licenses/LICENSE-2.0.html',
    url='https://github.com/radiasoft/rsfriction',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Environment :: Console',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: Apache Software License',
        'Natural Language :: English',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python',
        'Topic :: Utilities',
    ],
)

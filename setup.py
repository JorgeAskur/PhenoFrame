"""Mark wheels as platform-specific because phenoframe bundles a native library."""

from setuptools import Distribution, setup


class NativeLibraryDistribution(Distribution):
    def has_ext_modules(self) -> bool:
        return True


setup(distclass=NativeLibraryDistribution)

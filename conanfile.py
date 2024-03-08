import os
from conan import ConanFile
from conan.tools.cmake import CMakeToolchain, CMakeDeps


class MolppRecipe(ConanFile):
    settings = "os", "compiler", "build_type", "arch"

    def requirements(self):
        self.requires("eigen/3.4.0")
        self.requires("cpp-peglib/1.8.6")
        self.requires("gtest/1.14.0")

    def build_requirements(self):
        self.tool_requires("cmake/3.28.1")

    def generate(self):
        deps = CMakeDeps(self)
        deps.generate()
        tool_chain = CMakeToolchain(self)
        tool_chain.cache_variables["MOLPP_BUILD_TESTING"] = os.environ.get("BUILD_TESTING", "OFF")
        tool_chain.generate()

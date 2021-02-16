clear all
close all

Test = GraphOutput('Description','Test','Function',Type([sym('b1') sym('b2')],{sym('b1') sym('b2')},'b1*b2'),'Breakpoints',{GraphVertex_Internal.empty GraphVertex.empty});
% Test.Breakpoints = {GraphVertex_Internal.empty GraphVertex.empty} 